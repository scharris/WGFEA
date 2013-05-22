#module TMesh
#export TMesh

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

typealias Point Vec2
typealias PointNum Uint
typealias ElTag Uint32

# Element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELTYPE = 2
const TOKN_NUMTAGS = 3

# Gmsh-defined triangle type codes
const ELTYPE_3_NODE_TRIANGLE = 2
const ELTYPE_6_NODE_TRIANGLE = 9 # 3 vertex nodes and 3 on on edges

# Numbering the Sides of a Triangle
# In order to enumerate the sides of a triangle, we first provide a consistent way
# of ordering its vertex nodes. Vertex nodes are numbered by ordering the coordinate
# components of the nodes lexicographically in ascending order, ie. by first
# components most significantly then by second components for equal first components.
#
# Sides are then numbered as follows:
#     side 1: vertex 1 to vertex 2
#     side 2: vertex 1 to vertex 3
#     side 3: vertex 2 to vertex 3
#
#       O-v2                  O-v3
#   s1-   -s3          s2-    -s3
#  v1-O   O-v3              O-v2
#       ^s2        v1-O  ^s1
#
# sij Directed Side Notation
# Using the same vertex numbering, we will also use the notation sij to indicate the
# vector from vertex i to vertex j.

immutable RefTri
  el_type::Uint8
  s12::Vec2,
  s13::Vec2,
  dep_dims_by_sidenum::Array{Dim, 1}
  outward_normals_by_sidenum::Array{Vec2, 1}
end

immutable ElTri
  oshape::OShapeNum
  vertex_1::Point
end


type OShapeStore
  oshapenums_by_vertex_offsets::Dict{(Vec2,Vec2), OShapeNum}
  ref_tris_by_oshapenum::Array{RefTri,1}

  function OShapeStore(exp_num_oshapes::Int)
    const ref_tris = Array(RefTri,0)
    sizehint(ref_tris, exp_num_oshapes)

    const oshapenums_by_vertex_offsets = Dict{(Vec2,Vec2), OShapeNum}()
    sizehint(oshapenums_by_vertex_offsets, exp_num_oshapes)

    new(ref_tris, oshapenums_by_vertex_offsets)
  end
end


function find_or_create_oshape(pts::Array{Point,1}, el_tags::Array{AsciiString,1}, oshape_store::OShapeStore)
  const s12 = pts[2]-pts[1]
  const s13 = pts[3]-pts[1]
  # TODO: How to best do this?  Side vectors probably won't often be found exactly.
  #       Probably should mark as physical regions the initial mesh triangles
  #       within gmsh, via a new script which converts a .msh to .geo with the
  #       physical regions declared.  Then can progressively refine this .geo
  #       by splitting.  Each physical region should then only contain two kinds
  #       of triangles, one an inversion of the other.
  # Or, could just try looking for the point +- 0,1,2 eps in each coordinate of the
  # key pair, ie. do multiple hash lookups.
  const existing_oshapenum = get(oshapenums_by_vertex_offsets, (s12, s13), nothing)
  if existing_oshapenum != nothing
    existing_oshapenum::OShapeNum
  else
    const s23 = pts[3]-pts[2]
    const ref_tri = RefTri(el_type,
                           s12,
                           s13,
                           side_dep_dims(s12, s13, s23),
                           outward_normals(s12, s13, s23))
     # TODO
  end
end


immutable TriMesh

  fes::Array{ElTri,1}

  nbsideincls_by_nbsidenum::Array{NBSideInclusions,1}
  nbsidenums_by_feface::Dict{(FENum, FEFaceNum), NBSideNum}

  num_b_sides::Int64

  function TriMesh(ios::IOStream, est_ref_tris::Int)

    const pts_by_ptnum = read_points(ios)

    if !read_through_line(ios, "\$Elements\n")
      error("Could not find beginning of elements section in mesh file.")
    else
      # Estimate the number of finite elements.
      const decl_num_els = int(readline(ios)) # can include unwanted lower order elements in the count
      # Skip non-triangle elements which are assumed to be lower order, subtract skipped elements from declared number of elements from the file.
      line, lines_read = read_until(ios, is_triangle_el_or_endmarker)
      const est_num_fes = decl_num_els - (lines_read - 1)

      const fes = Array(ElTri, 0)
      sizehint(fes, est_num_fes)

      const side_fe_incls = Dict{(PointNum, PointNum), Array{(FENum, FEFaceNum),1}}()
      sizehint(side_fe_incls, est_num_fes * 3/2)

      const oshape_store = OShapeStore(est_ref_tris)

      num_nb_sides = 0
      num_b_sides = 0

      while line != "\$EndElements" && line != ""
        const toks = split(line, ' ')
        const el_type = int(toks[TOKN_ELTYPE])

        if el_type == ELTYPE_6_NODE_TRIANGLE
          error("6 node triangles are not yet supported.")

        elseif el_type == ELTYPE_3_NODE_TRIANGLE
          const el_pts = sort(map(ptnum_str -> pts_by_ptnum[int(ptnum_str)],
                                  toks[TOKN_NUMTAGS+int(toks[TOKN_NUMTAGS])+1:]))

          # Lookup or create the oriented shape for this element.
          const oshape = let el_tags = toks[TOKN_NUMTAGS+1 : TOKN_NUMTAGS+int(toks[TOKN_NUMTAGS])]
            find_or_create_oshape(el_pts, el_tags)
          end

          const el_tri = ElTri(oshape, el_pts[1])
          fes.push!(el_tri)
          const fenum = length(fes)

          # TODO
          # For each side, add to side_fe_incls an entry
          #   (lesser pt#, greater pt#) -> [(fe, rface)]
          #   if rhs has length 1 after adding, bump num_b_sides.
          #   if rhs has length 2 after adding, bump num_nb_sides, decrement num_b_sides.
        end

        line = readline(ios)
      end # line reading loop

      if line == ""
        error("Missing end of elements marker in mesh file.")
      end

      # TODO: Process side_incls, build nb sides data structures
      # nbsideincls_by_nbsidenum::Array{NBSideInclusions,1}
      # nbsidenums_by_feface::Dict{(FENum, FEFaceNum), NBSideNum}

      "TODO"
  end
  end
end

function side_dep_dims(side_1::Vec2, side_2::Vec2, side_3::Vec2)
  const dep_dims = Array(Dim, 0)
  for sv in (side_1, side_2, side_3)
    dep_dims.push(abs(sv[1]) > abs(sv[2]) ? dim(2) : dim(1))
  end
end

# Computing Outward Normals
# Let s be the z component of s12 x s13 (cross product).
# Then the outward normals are
#  side 1:  s12 x (0,0,s)
#  side 2: -s13 x (0,0,s)
#  side 3:  s23 x (0,0,s)
function outward_normals(s12::Vec2, s13::Vec2, s23::Vec2)
  # TODO
end


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
      const pt = (convert(R, float64(toks[2])), convert(R, float64(toks[3])))
      if length(toks) >= 4 && toks[4] != "0" && toks[4] != "0.0"
        error("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.")
      end
      pts[int(toks[1])] = pt
      l = strip(readline(ios))
    end
    if l == "" error("End of nodes section not found in mesh file.") end
    pts
  end
end

function read_through_line(io::IOStream, line::ASCIIString)
  l = readline(io)
  while l != "" && l != line
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted.
# Return either the line passing the condition, or "" if the stream was exhausted, and (in either case)
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

function is_triangle_el_or_endmarker(line::ASCIIString)
  if line == "\$EndElements"
    true
  else
    const toks = split(line, ' ')
    const el_type = int(toks[TOKN_ELTYPE])
    el_type == ELTYPE_3_NODE_TRIANGLE || el_type == ELTYPE_6_NODE_TRIANGLE
  end
end


immutable Vec2
  x::R,
  y::R
end

import Base.(+)
+(v1::Vec2, v2::Vec2) = Vec2(v1.x + v2.x, v1.y + v2.y)

import Base.(-)
-(v1::Vec2, v2::Vec2) = Vec2(v1.x - v2.x, v1.y - v2.y)

import Base.norm
norm(v::Vec2) = sqrt(v1.x*v1.x + v1.y*v1.y)

#end # end of module