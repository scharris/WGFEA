#module TMesh
#export TMesh

# require("TMesh"); ios = open("meshes/one_subdivided_triangle_mesh.msh"); tmsh = TriMesh(ios, 100)

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum, Mesh.fenum, Mesh.oshapenum, Mesh.nbsidenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

# vector type for representing relative node offsets and points in the mesh
immutable Vec
  _1::R
  _2::R
end
typealias Point Vec

# Gmsh element type code
typealias ElTypeNum Uint8
eltypenum(i::Integer) = convert(ElTypeNum, i)
eltypenum(s::String) = convert(ElTypeNum, int(s))

# reference triangle
immutable RefTri
  el_type::ElTypeNum
  v12::Vec
  v13::Vec
  dep_dims_by_sidenum::Array{Dim, 1}
  outward_normals_by_sidenum::Array{Vec, 1}
end

# finite element triangle
immutable ElTri
  oshape::OShapeNum
  vertex_1::Point
end

# The oriented shape store associates reference triangles with oriented shape numbers, and oriented shape
# numbers with triangle characteristics.  This enables an element's oriented shape and reference triangle
# to be retrieved when given its nodes.
type OShapeStore
  ref_tris_by_oshapenum::Array{RefTri,1}
  oshapenums_by_vertex_offsets::Dict{(Vec,Vec,ElTypeNum), OShapeNum}
  function OShapeStore(exp_num_oshapes::Integer)
    const ref_tris = sizehint(Array(RefTri,0), exp_num_oshapes)
    const oshapenums_by_vertex_offsets = sizehint(Dict{(Vec,Vec,ElTypeNum), OShapeNum}(), exp_num_oshapes)
    new(ref_tris, oshapenums_by_vertex_offsets)
  end
end


immutable TriMesh <: AbstractMesh

  # finite elements
  fes::Array{ElTri,1}

  # reference elements
  refs::Array{RefTri,1}

  # non-boundary side numbers by fe face.
  nbsidenums_by_feface::Dict{(FENum, FEFaceNum), NBSideNum}

  # inclusions of non-boundary sides in fe's, indexed by non-boundary side number
  nbsideincls_by_nbsidenum::Array{NBSideInclusions,1}

  # number of boundary sides
  num_b_sides::Int64

  # mesh construction from input stream
  function TriMesh(ios::IOStream, est_ref_tris::Int)
    const oshape_store = OShapeStore(est_ref_tris)

    const pts_by_nodenum = read_points(ios)

    if !read_through_line(ios, "\$Elements\n") error("Could not find beginning of elements section in mesh file.") end

    # Read to first triangle, estimating number of elements from that declared and number of lower order elements skipped.
    const decl_num_els = int(readline(ios)) # can include unwanted lower order elements in the count
    line, lines_read = read_until(ios, is_triangle_el_or_endmarker)
    const est_num_fes = decl_num_els - (lines_read - 1)
    println("Estimated $est_num_fes elements, $(lines_read - 1) lower order elements skipped.")

    const fes = sizehint(Array(ElTri, 0), est_num_fes)
    const side_endpoints_fe_incls = sizehint(Dict{(Point, Point), Array{(FENum, FEFaceNum),1}}(), uint64(est_num_fes * 3/2)) # 3+ sides per triangle, most included in 2 fes
    const verts = Array(Point, 3)
    num_nb_sides = 0
    num_b_sides = 0

    while line != "\$EndElements\n" && line != ""
      const toks = split(line, ' ')
      const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
      if is_lower_order_el_type(el_type) continue end
      if !is_triangle_el_type(el_type) error("Element type $el_type is not supported.") end
      const el_tags = toks[TOKN_ELLINE_NUMTAGS+1 : TOKN_ELLINE_NUMTAGS+int(toks[TOKN_ELLINE_NUMTAGS])]

      println("Filling vertexes for line $line.")

      # Get vertexes. Remaining points are not read because they are determined by the element type and vertexes.
      fill_verts_left_lower_first(toks, pts_by_nodenum, verts)

      println("Filled vertexes")

      const oshape = find_or_create_oshape(el_type, el_tags, verts, oshape_store)

      const el_tri = ElTri(oshape, verts[1])
      push!(fes, el_tri)
      const fe_num = fenum(length(fes))

      const nb_delta, b_delta = register_fe_faces_by_side_endpoints(fe_num, el_type, verts, side_endpoints_fe_incls)
      num_nb_sides += nb_delta
      num_b_sides += b_delta

      line = readline(ios)
    end # element line reading loop

    if line != "\$EndElements\n" error("Missing end of elements marker in mesh file.") end

    const nbsidenums_by_feface, nbsideincls_by_nbsidenum =
      create_nb_sides_data(side_endpoints_fe_incls, num_nb_sides_hint=num_nb_sides)

    new(fes,
        oshape_store.ref_tris_by_oshapenum,
        nbsidenums_by_feface,
        nbsideincls_by_nbsidenum,
        num_b_sides)
  end
end

# Fill the passed buffer with vertex points.  The vertexes are filled with the left-lower-most vertex first,
# otherwise matching (cyclically) the order in the mesh file.
function fill_verts_left_lower_first(toks::Array{String,1},
                                     pts_by_nodenum::Array{Point,1},
                                     verts_buf::Array{Point,1})
  const last_tag_tokn = TOKN_ELLINE_NUMTAGS + int(toks[TOKN_ELLINE_NUMTAGS])
  for i=1:3
    verts_buf[i] = pts_by_nodenum[uint64(toks[last_tag_tokn + i])]
  end

  # Rearrange the vertexes if necessary to make the first point the least lexicographically.
  if isless(verts_buf[2],verts_buf[1]) && isless(verts_buf[2],verts_buf[3])
    const pt_1 = verts_buf[1]
    verts_buf[1] = verts_buf[2]
    verts_buf[2] = verts_buf[3]
    verts_buf[3] = pt_1
  elseif isless(verts_buf[3],verts_buf[1]) && isless(verts_buf[3],verts_buf[2])
    const pt_3 = verts_buf[3]
    verts_buf[3] = verts_buf[2]
    verts_buf[2] = verts_buf[1]
    verts_buf[1] = pt_3
  end
end

# Register the finite element's side faces by their endpoints in the passed registry.
function register_fe_faces_by_side_endpoints(fe_num::FENum,
                                             el_type::ElTypeNum,
                                             verts::Array{Point,1},
                                             side_endpoints_fe_incls::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}})
  const side_endpoint_pairs = sideface_endpoint_pairs(el_type, verts, lesser_endpoints_first=true)
  nb_sides_delta = 0
  b_sides_delta = 0

  # Record the side faces included by this element.
  for sf=fefacenum(1):fefacenum(length(side_endpoint_pairs))
    const side_endpoints_stdord = side_endpoint_pairs[sf]
    const existing_side_incls = get(side_endpoints_fe_incls, side_endpoints_stdord, nothing)
    if existing_side_incls != nothing
      push!(existing_side_incls, (fe_num, sf))
      b_sides_delta -= 1
      nb_sides_delta += 1
    else
      side_endpoints_fe_incls[side_endpoints_stdord] = [(fe_num, sf)]
      b_sides_delta += 1
    end
  end
  nb_sides_delta, b_sides_delta
end

# Return the pair consisting of
#   1) nb side numbers by (fe,face) :: Dict{(FENum, FEFaceNum), NBSideNum},
#   2) an array of NBSideInclusions indexed by nb side number.
function create_nb_sides_data(side_endpoints_fe_incls::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}};
                              num_nb_sides_hint::Integer=0)
  const nbsidenums_by_feface = sizehint(Dict{(FENum, FEFaceNum), NBSideNum}(), 2*num_nb_sides_hint)
  const nbsideincls_by_nbsidenum = sizehint(Array(NBSideInclusions, 0), num_nb_sides_hint)

  for side_endpoints in keys(side_endpoints_fe_incls)
    const fe_faces = side_endpoints_fe_incls[side_endpoints]
    const num_fe_incls = length(fe_faces)
    if num_fe_incls == 2
      sort!(fe_faces)
      const (fe_1, side_face_1) = fe_faces[1]
      const (fe_2, side_face_2) = fe_faces[2]

      # Assign a new non-boundary side number.
      const nb_side_num = nbsidenum(length(nbsideincls_by_nbsidenum) + 1)

      # Add nb side inclusions structure for this nb side number.
      push!(nbsideincls_by_nbsidenum,
            NBSideInclusions(nb_side_num,
                             fe_1, side_face_1,
                             fe_2, side_face_2))

      # Map the two fe/face pairs to this nb side number.
      nbsidenums_by_feface[fe_faces[1]] = nb_side_num
      nbsidenums_by_feface[fe_faces[2]] = nb_side_num
    elseif num_fe_incls > 2
      error("Invalid mesh: side encountered with more than two including finite elements.")
    end
  end

  nbsidenums_by_feface, nbsideincls_by_nbsidenum
end

function side_dep_dims(el_type::ElTypeNum, verts::Array{Point,1})
  const sfs_btw_verts = num_side_faces_between_any_two_vertexes(el_type)
  const dep_dims = sizehint(Array(Dim, 0), 3*sfs_btw_verts)
  for vp=1:3  # vertex pairs
    const v1 = verts[vp]
    const v2 = verts[mod(vp,3)+1]
    const dep_dim = abs(v2._1-v1._1) > abs(v2._2-v1._2) ? dim(2) : dim(1)
    for sf=1:sfs_btw_verts
      push!(dep_dims, dep_dim)
    end
  end
  dep_dims
end

# Computing Outward Normals for Side Faces
# The outward normal for an inter-vertex vector w (extended to R^3 via 0 third component) is
#   n = (w x (0,0,cc)) / |w|
#     = cc (w[2], -w[1], 0) / |w|.
# where cc = 1 if w is part of a clockwise traversal of the triangle, and -1 otherwise.
# For any pair of successive inter-vertex vectors w_1 and w_2, cc can be computed as:
#   cc = sgn(z(w_1 x w_2)), where z is projection of component 3.
function outward_normals(el_type::ElTypeNum, verts::Array{Point,1})
  const sfs_btw_verts = num_side_faces_between_any_two_vertexes(el_type)
  const normals = sizehint(Array(Vec, 0), 3*sfs_btw_verts)
  # inter-vertex vectors
  const v12 = verts[2]-verts[1]
  const v23 = verts[3]-verts[2]
  const v31 = verts[1]-verts[3]
  const cc = counterclockwise(v12, v23) ? 1. : -1.
  for ivv in (v12, v23, v31)
    const len = norm(ivv)
    const n = Vec(cc * ivv._2/len, -cc * ivv._1/len)
    for sf=1:sfs_btw_verts
      push!(normals, n)
    end
  end
  normals
end


# oriented shape storage operations

function find_or_create_oshape(el_type::ElTypeNum,
                               el_tags::Array{String,1},
                               verts::Array{Point,1},
                               oshape_store::OShapeStore)
  # The reference oriented shape is determined by the vectors from node 1 to the other nodes, taken
  # as a pair in counterclockwise orientation, together with the mesh's element type code (which
  # determines how many additional nodes may be on each triangle side).
  const v12 = verts[2]-verts[1]
  const v13 = verts[3]-verts[1]
  const oshape_key = counterclockwise(v12, v13) ? (v12, v13, el_type) : (v13, v12, el_type)

  # TODO: This simple lookup is likely to often fail because of small floating point differences, so
  #       should try looking for the point +- 0,1,2 eps in each coordinate of the key pair, ie. do multiple hash lookups.
  const existing_oshape = get(oshape_store.oshapenums_by_vertex_offsets, oshape_key, nothing)

  if existing_oshape != nothing
    existing_oshape
  else
    const ref_tri = RefTri(el_type,
                           v12,
                           v13,
                           side_dep_dims(el_type, verts),
                           outward_normals(el_type, verts))

    push!(oshape_store.ref_tris_by_oshapenum, ref_tri)
    const oshape = oshapenum(length(oshape_store.ref_tris_by_oshapenum))
    oshape_store.oshapenums_by_vertex_offsets[oshape_key] = oshape
    oshape
  end
end

# element types


# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
const ELTYPE_6_NODE_TRIANGLE = eltypenum(9) # 3 vertex nodes and 3 on on edges
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)


function is_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE || el_type == ELTYPE_6_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

# This function defines the side faces enumeration for a finite element of given vertexes and mesh
# element type. Side faces are returned as an array of side endpoint pairs indexed by side face
# number. If lesser_endpoints_first is true, then each endpoint pair endpoint will have the lesser
# point (compared lexicographically) in the first component of the pair.
function sideface_endpoint_pairs(el_type::ElTypeNum, verts::Array{Point,1}; lesser_endpoints_first::Bool=false)
  const mk_pair = if lesser_endpoints_first; (pt1::Point, pt2::Point) -> isless(pt1,pt2) ? (pt1, pt2) : (pt2, pt1)
                  else (pt1::Point, pt2::Point) -> (pt1, pt2) end
  if el_type == ELTYPE_3_NODE_TRIANGLE
    [mk_pair(verts[1],verts[2]), mk_pair(verts[2],verts[3]), mk_pair(verts[3],verts[1])]
  elseif el_type == ELTYPE_6_NODE_TRIANGLE
    error("6 node triangles not yet supported.")
  else
    error("Unsupported element type $el_type.")
  end
end

pair_lesser_fst(pt1::Point, pt2::Point) = if isless(pt1,pt2) (pt1, pt2) else (pt2, pt1) end

function num_side_faces_between_any_two_vertexes(el_type::ElTypeNum)
  if el_type == ELTYPE_3_NODE_TRIANGLE
    1
  elseif el_type == ELTYPE_6_NODE_TRIANGLE
    2
  else
    error("Unsupported element type $el_type.")
  end
end

# file reading utilities

function read_points(ios::IOStream)
  if !read_through_line(ios, "\$Nodes\n")
    error("Nodes section not found in mesh input file.")
  else
    # Next line should be the vert count.
    const count = int(readline(ios))
    const pts = Array(Point, count)
    l = strip(readline(ios))
    while l != "\$EndNodes" && l != ""
      const toks = split(l, ' ')
      const pt = Point(convert(R, float64(toks[TOKN_NODELINE_POINT1])), convert(R, float64(toks[TOKN_NODELINE_POINT2])))
      if length(toks) >= TOKN_NODELINE_POINT3 && toks[TOKN_NODELINE_POINT3] != "0" && toks[TOKN_NODELINE_POINT3] != "0.0"
        error("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.")
      end
      pts[uint64(toks[TOKN_NODELINE_ELNUM])] = pt
      l = strip(readline(ios))
    end
    if l == "" error("End of nodes section not found in mesh file.") end
    pts
  end
end

function is_triangle_el_or_endmarker(line::ASCIIString)
  if line == "\$EndElements"
    true
  else
    const toks = split(line, ' ')
    is_triangle_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
  end
end

function read_through_line(io::IOStream, line::ASCIIString)
  l = readline(io)
  while l != "" && l != line
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IOStream, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end


# vector operations

import Base.getindex
getindex(v::Vec, i::Integer) = if i == 1 v._1 elseif i == 2 v._2 else throw(BoundsError()) end

import Base.(+)
+(v::Vec, w::Vec) = Vec(v._1 + w._1, v._2 + w._2)

import Base.(-)
-(v::Vec, w::Vec) = Vec(v._1 - w._1, v._2 - w._2)

import Base.norm
norm(v::Vec) = hypot(v._1, v._2)

import Base.isless
isless(v::Vec, w::Vec) = v._1 < w._1 || v._1 == w._1 && v._2 < w._2

counterclockwise(u::Vec, v::Vec) =
  u._1*v._2 - u._2*v._1 > 0


# constants related to the Gmsh format

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3

#end # end of module