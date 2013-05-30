module TMesh
export TriMesh, Vec, RefTri, ElTri

# using TMesh; ios = open("meshes/one_subdivided_triangle_mesh.msh"); tmsh = TriMesh(ios, 100)

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum, Mesh.fenum, Mesh.nbsidenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

# vector type for representing relative node offsets and points in the mesh
immutable Vec
  _1::R
  _2::R
end
typealias Point Vec

# reference triangle
immutable RefTri
  v12::Vec
  v13::Vec
  nums_faces_between_vertexes::(Int,Int,Int) # v1v2, v2v3, v3v1
  outward_normals_by_sideface::Array{Vec,1}
  diameter_inv::R
end

# finite element triangle
immutable ElTri
  oshape::OShapeNum
  v1::Point
end

immutable TriMesh <: AbstractMesh
  # finite elements, indexed by finite element number
  fes::Array{ElTri,1}

  # reference triangles, indexed by oriented shape number
  oshapes::Array{RefTri,1}

  # non-boundary side numbers by fe face.
  nbsidenums_by_feface::Dict{(FENum,FEFaceNum), NBSideNum}

  # inclusions of non-boundary sides in fe's, indexed by non-boundary side number
  nbsideincls_by_nbsidenum::Array{NBSideInclusions,1}

  # dependent dimensions for nb sides
  # For a given non-boundary side, this is a coordinate that is a function of the others over the side.
  dep_dims_by_nbsidenum::Array{Dim,1}

  # number of boundary sides
  num_b_sides::Int64
end


# Construct from Gmsh .msh formatted input stream, and number of subdivision operations to perform within each
# mesh element from the stream.  Some elements may have additional iterations specified via Gmsh tags.
function TriMesh(ios::IO, base_subdiv_iters::Integer)
  const base_subdiv_iters = base_subdiv_iters >= 0 ? uint(base_subdiv_iters) : throw(ArgumentError("Number of iterations must be non-negative."))

  # Read the points from the input mesh, storing by mesh node number.  The elements section of the input
  # stream will refer to these node numbers.
  const mesh_pts_by_nodenum = read_points(ios)

  # Read to the beginning of the elements section.
  if !read_through_line(ios, "\$Elements\n") error("Could not find beginning of elements section in mesh file.") end

  # Skip point and line elements, estimating number of input mesh triangles from the number of elements declared.
  const decl_num_els = uint64(readline(ios)) # This count can include unwanted lower order elements.
  line, lines_read = read_until(ios, is_polytope_or_endmarker)
  const est_mesh_tris = decl_num_els - (lines_read - 1)

  # Make estimates of the number of finite elements and reference triangles for storage allocation.
  # This estimate ignores optional additional subdivisions that may be specified for some input elements.
  const est_fes = est_mesh_tris * 4^base_subdiv_iters
  # We'll estimate 2 reference triangles for each mesh element if we are subdividing, 1 if not. We ignore
  # for this estimate any additional reference triangles for hanging node support.
  const est_ref_tris = (base_subdiv_iters > 0 ? 2 : 1) * est_mesh_tris
  println("Estimating $(int(est_fes)) finite elements, $(int(est_ref_tris)) reference triangles.")

  # Data to be updated as input element lines are processed.
  const oshapes = sizehint(Array(RefTri,0), est_ref_tris)
  const fes = sizehint(Array(ElTri, 0), est_fes)
  const fe_faces_by_endpoints = sizehint(Dict{(Point,Point), Array{(FENum,FEFaceNum),1}}(), uint64(est_fes * 3/2)) # estimate assumes most fes have 3 sides, each in 2 fes
  num_nb_sides = 0
  num_b_sides = 0

  # Function via which all finite elements will be created during the processing of input element lines.
  const fe_maker = function(oshapenum::OShapeNum, v1::Point, v2::Point, v3::Point)
     push!(fes, ElTri(oshapenum, v1))
     const fe_num = fenum(length(fes))
     const nums_faces_btw_verts = oshapes[oshapenum].nums_faces_between_vertexes
     const nb_delta, b_delta = register_fe_faces_by_endpoints(fe_num, v1,v2,v3, nums_faces_btw_verts, fe_faces_by_endpoints)
     num_nb_sides += nb_delta
     num_b_sides += b_delta
  end

  # process input mesh elements
  while line != "\$EndElements\n" && line != ""
    process_el_line(line, mesh_pts_by_nodenum, base_subdiv_iters, fe_maker, oshapes)
    line = readline(ios)
  end

  if line != "\$EndElements\n" error("Missing end of elements marker in mesh file.") end

  # Create the final non-boundary sides data structures based on the mapping of side endpoints to fe faces.
  const nbsidenums_by_feface, nbsideincls_by_nbsidenum, dep_dims_by_nbsidenum =
    create_nb_sides_data(fe_faces_by_endpoints, num_nb_sides_hint=num_nb_sides)

  TriMesh(fes,
          oshapes,
          nbsidenums_by_feface,
          nbsideincls_by_nbsidenum,
          dep_dims_by_nbsidenum,
          num_b_sides)
end


function process_el_line(line::String,
                         mesh_pts_by_nodenum::Array{Point,1},
                         base_subdiv_iters::Uint,
                         fe_maker::Function,
                         oshapes::Array{RefTri,1}) # will be appended to
  const toks = split(line, ' ')

  const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
  if is_lower_order_el_type(el_type) return end
  if !is_3_node_triangle_el_type(el_type) error("Element type $el_type is not supported.") end

  const el_tags = toks[TOKN_ELLINE_NUMTAGS+1 : TOKN_ELLINE_NUMTAGS+int(toks[TOKN_ELLINE_NUMTAGS])]

  const v1,v2,v3 = lookup_vertexes(toks, mesh_pts_by_nodenum)

  # Add any extra subdivision iterations to be done in this element.
  const subdiv_iters = base_subdiv_iters + extra_subdiv_iters(v1,v2,v3, el_tags)

  # For each of the three sides of this mesh element to be subdivided, the generated subdivision
  # elements can be made to have more than one face between vertexes lying on the side.  This is
  # used to support "hanging" nodes where a finer subdivision is adjacent to this one.  The triplet
  # returned represents the number of faces between element vertex pairs embedded in the input mesh
  # element sides from v1 to v2, v2 to v3, and v3 to v1, respectively.
  const nums_faces_btw_verts = nums_faces_between_vertexes(v1,v2,v3, el_tags)

  # Register the primary reference triangles for our mesh element's subdivisions.
  const pri_oshapenums_by_nums_faces_btw_verts = register_primary_ref_tris(v1,v2,v3, nums_faces_btw_verts, subdiv_iters, oshapes)

  if subdiv_iters > 0
    # Register the secondary reference triangle.
    push!(oshapes, secondary_ref_tri(v1,v2,v3, subdiv_iters))
    const sec_oshapenum = Mesh.oshapenum(length(oshapes))

    # Do the subdivisions.
    subdivide_primary(v1,v2,v3,
                      subdiv_iters,
                      nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)
  else # no subdivision to be done
    # The mesh element itself is our finite element, with its own reference triangle (added in register_primary_ref_tris above).
    const oshapenum = pri_oshapenums_by_nums_faces_btw_verts(nums_faces_btw_verts)
    fe_maker(oshapenum, v1,v2,v3)
  end
end # process_el_line


# Create primary reference triangles for the given triangle to be subdivided. These reference triangles
# are rescaled translations of the original undivided triangle. If nums_faces_btw_verts is other than
# (1,1,1), then multiple reference triangles will be generated, differing in the number of side faces
# between their vertexes for support of "hanging" nodes.  The function returns a function mapping the
# numbers of side faces between vertexes to the newly registered reference triangle (oshape) number.
function register_primary_ref_tris(v1::Point, v2::Point, v3::Point,
                                   nums_faces_btw_verts::(Int,Int,Int),
                                   subdiv_iters::Integer,
                                   oshapes::Array{RefTri,1})
  num_added = 0
  for nsf1 in subdiv_iters == 0 || nums_faces_btw_verts[1] == 1 ? nums_faces_btw_verts[1] : (1, nums_faces_btw_verts[1]),
      nsf2 in subdiv_iters == 0 || nums_faces_btw_verts[2] == 1 ? nums_faces_btw_verts[2] : (1, nums_faces_btw_verts[2]),
      nsf3 in subdiv_iters == 0 || nums_faces_btw_verts[3] == 1 ? nums_faces_btw_verts[3] : (1, nums_faces_btw_verts[3])
    push!(oshapes, primary_ref_tri(v1,v2,v3, (nsf1,nsf2,nsf3), subdiv_iters))
    num_added += 1
  end

  # Return a lookup function for these new primary oshape numbers by numbers of side faces between vertexes.
  const first_new_oshapenum = length(oshapes)-(num_added-1)
  const last_new_oshapenum = length(oshapes)
  function(nums_faces_btw_verts::(Int,Int,Int))
    for os=first_new_oshapenum:last_new_oshapenum
      if oshapes[os].nums_faces_between_vertexes == nums_faces_btw_verts
        return Mesh.oshapenum(os)
      end
    end
    error("Reference triangle not found by numbers of side faces between vertexes: $nums_faces_btw_verts.")
  end
end

function primary_ref_tri(v1::Point, v2::Point, v3::Point,
                         nums_faces_btw_verts::(Int,Int,Int),
                         subdiv_iters::Integer)
  const scale = 1./2^subdiv_iters
  const scaled_v12 = scale*(v2-v1)
  const scaled_v13 = scale*(v3-v1)
  const norm_scaled_v23 = scale * norm(v3-v2)
  const diameter_inv = 1./max(norm(scaled_v12), norm(scaled_v13), norm_scaled_v23)
  RefTri(scaled_v12,
         scaled_v13,
         nums_faces_btw_verts,
         outward_normals(v1,v2,v3, nums_faces_btw_verts),
         diameter_inv)
end

# Returns the reference element for the secondary, less numerous, inverted triangles which
# form in the middle positions after triangle subdivisions. Note that passed vertexes represent
# those of an outermost triangle to be subdivided, and hence describe a triangle of the primary,
# not secondary, oriented shape.
function secondary_ref_tri(pri_v1::Point, pri_v2::Point, pri_v3::Point, # vertexes of the undivided, *primary* shape
                           subdiv_iters::Integer)
  if subdiv_iters == 0 error("Secondary reference triangle requires non-zero number of iterations.") end
  const v1,v2,v3 = .5(pri_v1+pri_v2), .5(pri_v2+pri_v3), .5(pri_v3+pri_v1)
  const scale = 1./2^(subdiv_iters-1) # v1,v2,v3 already represents one subdivision
  const scaled_v12 = scale*(v2-v1)
  const scaled_v13 = scale*(v3-v1)
  const norm_scaled_v23 = scale * norm(v3-v2)
  const diameter_inv = 1./max(norm(scaled_v12), norm(scaled_v13), norm_scaled_v23)
  RefTri(scaled_v12,
         scaled_v13,
         ONE_FACE_BETWEEN_VERTEX_PAIRS,
         outward_normals(v1,v2,v3),
         diameter_inv)
end


function subdivide_primary(v1::Point, v2::Point, v3::Point,
                           iters::Uint,
                           nums_faces_btw_verts::(Int,Int,Int),
                           pri_oshapenums_by_nums_faces_btw_verts::Function, sec_oshapenum::OShapeNum,
                           fe_maker::Function)
  if iters == 0
    const oshapenum = pri_oshapenums_by_nums_faces_btw_verts(nums_faces_btw_verts)
    fe_maker(oshapenum, v1,v2,v3)
  else
    const midpt_12 = 0.5(v1+v2)
    const midpt_13 = 0.5(v1+v3)
    const midpt_23 = 0.5(v2+v3)

    # The sub-triangle including v1 has its first and last vertex pairs embedded in the original triangle's first and
    # third sides, and so inherit the corresponding numbers of faces between vertexes from nums_faces_btw_verts.
    # The middle vertex pair (midpt_12 to midpt_13) of this sub-triangle does not lie along an original side, so has
    # only one face between these vertexes.  Similar arguments apply for the remaining sub-triangles.
    const st1_nums_faces_btw_verts = (nums_faces_btw_verts[1], 1, nums_faces_btw_verts[3])
    subdivide_primary(v1, midpt_12, midpt_13,
                      iters-1,
                      st1_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    const st2_nums_faces_btw_verts = (nums_faces_btw_verts[1], nums_faces_btw_verts[2], 1)
    subdivide_primary(midpt_12, v2, midpt_23,
                      iters-1,
                      st2_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    const st3_nums_faces_btw_verts = (1, nums_faces_btw_verts[2], nums_faces_btw_verts[3])
    subdivide_primary(midpt_13, midpt_23, v3,
                      iters-1,
                      st3_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    # secondary sub-triangle
    subdivide_secondary(midpt_12, midpt_23, midpt_13,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)
  end
end

const ONE_FACE_BETWEEN_VERTEX_PAIRS = (1,1,1)

function subdivide_secondary(v1::Point, v2::Point, v3::Point,
                             iters::Uint,
                             pri_oshapenums_by_nums_faces_btw_verts::Function, sec_oshapenum::OShapeNum,
                             fe_maker::Function)
  if iters == 0
    fe_maker(sec_oshapenum, v1,v2,v3)
  else
    const midpt_12 = 0.5(v1+v2)
    const midpt_13 = 0.5(v1+v3)
    const midpt_23 = 0.5(v2+v3)

    subdivide_secondary(v1, midpt_12, midpt_13,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    subdivide_secondary(midpt_12, v2, midpt_23,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    subdivide_secondary(midpt_13, midpt_23, v3,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    # secondary sub-triangle
    subdivide_primary(midpt_13, midpt_12, midpt_23,
                      iters-1,
                      ONE_FACE_BETWEEN_VERTEX_PAIRS,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)
  end
end


function lookup_vertexes(toks::Array{String,1}, mesh_pts_by_nodenum::Array{Point,1})
  const last_tag_tokn = TOKN_ELLINE_NUMTAGS + int(toks[TOKN_ELLINE_NUMTAGS])
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 1])],
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 2])],
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 3])]
end


# Register the finite element's side faces by their endpoints in the passed registry.
function register_fe_faces_by_endpoints(fe_num::FENum,
                                        v1::Point, v2::Point, v3::Point,
                                        nums_faces_btw_verts::(Int,Int,Int),
                                        fe_faces_by_endpoints::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}})
  const sf_endpt_pairs = sideface_endpoint_pairs(v1,v2,v3, nums_faces_btw_verts, lesser_endpts_first=true)
  nb_sides_delta = 0
  b_sides_delta = 0

  # Record the side faces included by this elementsment.
  for sf=fefacenum(1):fefacenum(length(sf_endpt_pairs))
    const sf_endpts = sf_endpt_pairs[sf]
    const existing_side_incls = get(fe_faces_by_endpoints, sf_endpts, nothing)
    if existing_side_incls != nothing
      push!(existing_side_incls, (fe_num, sf))
      b_sides_delta -= 1
      nb_sides_delta += 1
    else
      fe_faces_by_endpoints[sf_endpts] = [(fe_num, sf)]
      b_sides_delta += 1
    end
  end
  nb_sides_delta, b_sides_delta
end

# Return the triplet consisting of
#   1) nb side numbers by (fe,face) :: Dict{(FENum, FEFaceNum), NBSideNum},
#   2) an array of NBSideInclusions indexed by nb side number.
#   3) an array of side dependent dimension by nb side number
function create_nb_sides_data(fe_faces_by_endpoints::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}};
                              num_nb_sides_hint::Integer=0)
  const nbsidenums_by_feface = sizehint(Dict{(FENum, FEFaceNum), NBSideNum}(), 2*num_nb_sides_hint)
  const nbsideincls_by_nbsidenum = sizehint(Array(NBSideInclusions, 0), num_nb_sides_hint)
  const dep_dims_by_nbsidenum = sizehint(Array(Dim, 0), num_nb_sides_hint)

  for side_endpoints in keys(fe_faces_by_endpoints)
    const fe_faces = fe_faces_by_endpoints[side_endpoints]
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

      const dep_dim = let
                        const ep1,ep2 = side_endpoints
                        abs(ep2._1-ep1._1) >= abs(ep2._2-ep1._2) ? dim(2) : dim(1)
                      end
      push!(dep_dims_by_nbsidenum, dep_dim)
    elseif num_fe_incls > 2
      error("Invalid mesh: side encountered with more than two including finite elements.")
    end
  end

  nbsidenums_by_feface, nbsideincls_by_nbsidenum, dep_dims_by_nbsidenum
end

# Computing Outward Normals for Side Faces
# The outward normal for an inter-vertex vector w (extended to R^3 via 0 third component) is
#   n = (w x (0,0,cc)) / |w|
#     = cc (w[2], -w[1], 0) / |w|.
# where cc = 1 if w is part of a clockwise traversal of the triangle, and -1 otherwise.
# For any pair of successive inter-vertex vectors w_1 and w_2, cc can be computed as:
#   cc = sgn(z(w_1 x w_2)), where z is projection of component 3.
function outward_normals(v1::Point, v2::Point, v3::Point, nums_faces_btw_verts::(Int,Int,Int)=(1,1,1))
  const normals = sizehint(Array(Vec, 0), sum(nums_faces_btw_verts))
  # inter-vertex vectors
  const v12 = v2-v1
  const v23 = v3-v2
  const v31 = v1-v3
  const cc = counterclockwise(v12, v23) ? 1. : -1.
  for (ivv, num_faces_btw_verts) in ((v12, nums_faces_btw_verts[1]),
                                     (v23, nums_faces_btw_verts[2]),
                                     (v31, nums_faces_btw_verts[3]))
    const len = norm(ivv)
    const n = Vec(cc * ivv._2/len, -cc * ivv._1/len)
    for sf=1:num_faces_btw_verts
      push!(normals, n)
    end
  end
  normals
end


# This function defines the side faces enumeration for a finite element of given vertexes and numbers
# of faces between vertexes. Side faces are returned as an array of side endpoint pairs indexed by
# side face number. If lesser_endpts_first is true, then each endpoint pair endpoint will have
# the lesser point (compared lexicographically) in the first component of the pair.
function sideface_endpoint_pairs(v1::Point, v2::Point, v3::Point,
                                 nums_faces_btw_verts::(Int,Int,Int);
                                 lesser_endpts_first::Bool=false)
  const num_side_faces = sum(nums_faces_btw_verts)
  const sf_endpt_pairs = sizehint(Array((Point,Point),0), num_side_faces)
  const vert_pairs = ((v1,v2), (v2,v3), (v3,v1))
  for i=1:3
    const va,vb = vert_pairs[i]
    const num_faces_btw_va_vb = nums_faces_btw_verts[i]
    if num_faces_btw_va_vb == 1
      push!(sf_endpt_pairs, mk_endpoint_pair(va, vb, lesser_endpts_first))
    elseif num_faces_btw_va_vb == 2
      const midpt = 0.5(va+vb)
      push!(sf_endpt_pairs, mk_endpoint_pair(va, midpt, lesser_endpts_first))
      push!(sf_endpt_pairs, mk_endpoint_pair(midpt, vb, lesser_endpts_first))
    else
      error("Only 1 or 2 faces between triangle vertexes are currently supported.")
    end
  end
  sf_endpt_pairs
end

function mk_endpoint_pair(pt1::Point, pt2::Point, lesser_pt_first::Bool)
  if lesser_pt_first isless(pt1, pt2) ? (pt1, pt2) : (pt2, pt1)
  else (pt1, pt2) end
end


function extra_subdiv_iters(v1::Point, v2::Point, v3::Point, el_tags::Array{String,1})
  # TODO: Allow specifying extra iterations by tagging a mesh element.
  0
end

# For each of the three sides of the indicated mesh element to be subdivided, we can specify here the
# number of faces that the generated subtriangles should have between their vertexes which lie on this side.
# This allows these subdivision elements to meet those of a finer subdivision in a mesh element adjacent
# to this element ("hanging nodes").  The numbers are returned as a triplet of integers, corresponding to
# the face counts to be generated for elements' sides within sides v1v2, v2v3, and v3v1.
function nums_faces_between_vertexes(v1::Point, v2::Point, v3::Point, el_tags::Array{String,1})
  ONE_FACE_BETWEEN_VERTEX_PAIRS
end

# file reading utilities

function read_points(ios::IO)
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

function is_polytope_or_endmarker(line::ASCIIString)
  if line == "\$EndElements"
    true
  else
    const toks = split(line, ' ')
    !is_lower_order_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
  end
end

function read_through_line(io::IO, line::ASCIIString)
  l = readline(io)
  while l != "" && l != line
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IO, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end


# vector operations

#import Base.getindex
#getindex(v::Vec, i::Integer) = if i == 1 v._1 elseif i == 2 v._2 else throw(BoundsError()) end

import Base.(+)
+(v::Vec, w::Vec) = Vec(v._1 + w._1, v._2 + w._2)

import Base.(-)
-(v::Vec, w::Vec) = Vec(v._1 - w._1, v._2 - w._2)

import Base.(*)
*(r::Real, v::Vec) = Vec(r*v._1, r*v._2)

import Base.norm
norm(v::Vec) = hypot(v._1, v._2)

import Base.isless
isless(v::Vec, w::Vec) = v._1 < w._1 || v._1 == w._1 && v._2 < w._2

import Base.dot
dot(v::Vec, w::Vec) = v._1 * w._1 + v._2 * w._2

counterclockwise(u::Vec, v::Vec) =
  u._1*v._2 - u._2*v._1 > 0


# definitions related to the Gmsh format

typealias ElTypeNum Uint8
eltypenum(i::Integer) = if i>=1 && i<=typemax(ElTypeNum) convert(ElTypeNum, i) else error("invalid element type: $i") end
eltypenum(s::String) = eltypenum(int(s))

# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)

function is_3_node_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3


# String Representations
import Base.string, Base.show, Base.print

# string representation for monomials
function string(m::TriMesh)
  "TriMesh: $(length(m.fes)) elements, $(length(m.oshapes)) reference elements, $(length(m.nbsideincls_by_nbsidenum)) non-boundary sides, $(m.num_b_sides) boundary sides"
end
print(io::IO, m::TriMesh) = print(io, string(m))
show(io::IO, m::TriMesh) = print(io, m)

end # end of module
