module RMesh
export RectMesh,
       MeshCoord, mesh_coord, mesh_ldim, mesh_ldims,
       lesser_side_face_perp_to_axis, greater_side_face_perp_to_axis,
       exportAsGmshSurface

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

# type of a single logical mesh coordinates component (column or row or stack, etc)
typealias MeshCoord Uint64
mesh_coord(i::Integer) = convert(MeshCoord, i)
mesh_ldim(i::Integer) = convert(MeshCoord, i)
mesh_ldims(dims::Integer...) = [convert(MeshCoord,l) for l in dims]


immutable NBSideGeom
  perp_axis::Dim
  mesh_coords::Array{MeshCoord, 1}
end

const default_integration_rel_err = 1e-12
const default_integration_abs_err = 1e-12

immutable RectMesh <: AbstractMesh

  space_dim::Dim

  # Mesh coordinate ranges in R^d defining the boundaries of the mesh.
  min_bounds::Array{R,1}
  max_bounds::Array{R,1}

  # Logical dimensions of the mesh in discrete mesh axis coordinates,
  # with directions corresponding to the coordinate axes (cols, rows,...).
  mesh_ldims::Array{MeshCoord,1}

  # actual dimensions of any single finite element
  fe_dims::Array{R,1}
  fe_dims_wo_dim::Array{Array{R,1},1}

  cumprods_mesh_ldims::Array{FENum,1}

  cumprods_nb_side_mesh_ldims_by_perp_axis::Array{Array{NBSideNum,1},1}

  first_nb_side_nums_by_perp_axis::Array{NBSideNum,1}

  num_fes::FENum
  num_nb_sides::NBSideNum
  num_side_faces_per_fe::FEFaceNum

  rect_diameter::R
  rect_diameter_inv::R

  one_mon::Monomial


  # integration support members
  intgd_args_work_array::Vector{R}
  space_dim_zeros::Vector{R}
  space_dim_less_one_zeros::Vector{R}
  integration_rel_err::R
  integration_abs_err::R

  function RectMesh(min_bounds::Array{R,1},
                    max_bounds::Array{R,1},
                    mesh_ldims::Array{MeshCoord,1},
                    integration_rel_err::R,
                    integration_abs_err::R)
    const space_dim = length(min_bounds) | uint
    assert(length(max_bounds) == space_dim, "min and max bound lengths should match")
    assert(length(mesh_ldims) == space_dim, "logical dimensions length does not match physical bounds length")

    const fe_dims = make_fe_dims(min_bounds, max_bounds, mesh_ldims)
    const fe_dims_wo_dim = make_fe_dims_with_drops(fe_dims)

    const cumprods_mesh_ldims = cumprod(mesh_ldims)
    const cumprods_nb_side_mesh_ldims = make_cumprods_nb_side_mesh_ldims_by_perp_axis(mesh_ldims)

    const nb_side_counts_by_perp_axis = map(last, cumprods_nb_side_mesh_ldims)
    const first_nb_side_nums_by_perp_axis = cumsum(vcat(1, nb_side_counts_by_perp_axis[1:space_dim-1]))

    const num_fes = last(cumprods_mesh_ldims)
    const num_nb_sides = sum(nb_side_counts_by_perp_axis)
    const num_side_faces_per_fe = fefacenum(2 * space_dim)

    const rect_diameter = sqrt(dot(fe_dims, fe_dims))
    const rect_diameter_inv = 1./rect_diameter

    new(dim(space_dim),
        min_bounds,
        max_bounds,
        mesh_ldims,
        fe_dims,
        fe_dims_wo_dim,
        cumprods_mesh_ldims,
        cumprods_nb_side_mesh_ldims,
        first_nb_side_nums_by_perp_axis,
        num_fes,
        num_nb_sides,
        num_side_faces_per_fe,
        rect_diameter,
        rect_diameter_inv,
        Monomial(zeros(Deg,space_dim)),
        Array(R, space_dim), # integrand args work array
        zeros(R, space_dim),
        zeros(R, space_dim-1),
        integration_rel_err,
        integration_abs_err)
  end
end # type RectMesh

RectMesh(min_bounds::Array{R,1},
         max_bounds::Array{R,1},
         mesh_ldims::Array{MeshCoord,1}) =
  RectMesh(min_bounds, max_bounds, mesh_ldims, default_integration_rel_err, default_integration_abs_err)

import Base.isequal
function isequal(mesh1::RectMesh, mesh2::RectMesh)
  if is(mesh1, mesh2)
    true
  else
    mesh1.min_bounds == mesh2.min_bounds &&
    mesh1.max_bounds == mesh2.max_bounds &&
    mesh1.mesh_ldims == mesh2.mesh_ldims
  end
end

import Base.hash
hash(mesh::RectMesh) = hash(mesh.min_bounds) + 3 * hash(mesh.max_bounds) + 5 * hash(mesh.mesh_ldims)

# Auxiliary construction functions

function make_fe_dims(min_bounds::Array{R,1}, max_bounds::Array{R,1}, mesh_ldims::Array{MeshCoord,1})
  const space_dim = length(min_bounds)
  const dims = Array(R, space_dim)
  for i=1:space_dim
    const bounds_diff = max_bounds[i] - min_bounds[i]
    const ldim_i = mesh_ldims[i]
    assert(bounds_diff > zeroR, "improper mesh bounds")
    assert(ldim_i > 0, "non-positive logical mesh dimension")
    dims[i] = bounds_diff/ldim_i
  end
  dims
end

function make_fe_dims_with_drops(fe_dims::Array{R,1})
  const d = length(fe_dims)
  const dims_wo_dim = Array(Array{R,1}, d)
  for r=1:d
    dims_wo_dim[r] = fe_dims[[1:r-1,r+1:d]]
  end
  dims_wo_dim
end

function make_cumprods_nb_side_mesh_ldims_by_perp_axis(fe_mesh_ldims::Array{MeshCoord,1})
  const space_dim = length(fe_mesh_ldims)
  [ [make_cumprod_nb_side_mesh_ldims_to(r, perp_axis, fe_mesh_ldims) for r=1:space_dim]
    for perp_axis=1:space_dim ]
end

function make_cumprod_nb_side_mesh_ldims_to(r::Int, perp_axis::Int, fe_mesh_ldims::Array{MeshCoord,1})
  prod = one(NBSideNum)
  for i=1:r
    prod *= (i != perp_axis ? fe_mesh_ldims[i] : fe_mesh_ldims[i]-1)
  end
  prod::NBSideNum
end


##############################################
## Implement functions required of all meshes.

import Mesh.space_dim
space_dim(mesh::RectMesh) = mesh.space_dim

import Mesh.one_mon
one_mon(mesh::RectMesh) = mesh.one_mon

import Mesh.num_fes
num_fes(mesh::RectMesh) = mesh.num_fes

import Mesh.num_nb_sides
num_nb_sides(mesh::RectMesh) = mesh.num_nb_sides

import Mesh.num_oriented_element_shapes
num_oriented_element_shapes(mesh::RectMesh) = Mesh.oshape_one

import Mesh.oriented_shape_for_fe
oriented_shape_for_fe(fe::FENum, mesh::RectMesh) = Mesh.oshape_one

import Mesh.num_side_faces_for_fe
num_side_faces_for_fe(fe::FENum, mesh::RectMesh) = mesh.num_side_faces_per_fe

import Mesh.num_side_faces_for_shape
num_side_faces_for_shape(oshape::OShapeNum, mesh::RectMesh) = mesh.num_side_faces_per_fe

#import Mesh.dependent_dim_for_nb_side
#dependent_dim_for_nb_side(i::NBSideNum, mesh::RectMesh) = perp_axis_for_nb_side(i, mesh)

import Mesh.dependent_dim_for_oshape_side
dependent_dim_for_oshape_side(fe_oshape::OShapeNum, side_face::FEFaceNum, mesh::RectMesh) =
  side_face_perp_axis(side_face)

import Mesh.fe_inclusions_of_nb_side
function fe_inclusions_of_nb_side(n::NBSideNum, mesh::RectMesh)
  const sgeom = nb_side_geom(n, mesh)
  const a = sgeom.perp_axis
  const lesser_fe = fe_with_mesh_coords(sgeom.mesh_coords, mesh)
  const greater_fe =  lesser_fe + (a == 1 ? 1 : mesh.cumprods_mesh_ldims[a-1])
  NBSideInclusions(n,
                   lesser_fe, greater_side_face_perp_to_axis(a),
                   greater_fe, lesser_side_face_perp_to_axis(a))
end

import Mesh.nb_side_num_for_fe_side
function nb_side_num_for_fe_side(fe::FENum, side_face::FEFaceNum, mesh::RectMesh)
  const a = side_face_perp_axis(side_face)
  const side_mesh_coords =
    if side_face_is_lesser_on_perp_axis(side_face)
      # Side is the lesser along perp axis, ie. fe is the greater one including the side: perp axis dim differs.
      let coords = fe_mesh_coords(fe, mesh)
        coords[a] -= 1
        coords
      end
    else
      fe_mesh_coords(fe, mesh)
    end
  nb_side_with_mesh_coords(side_mesh_coords, a, mesh)
end


import Mesh.is_boundary_side
function is_boundary_side(fe::FENum, side_face::FEFaceNum, mesh::RectMesh)
  const a = side_face_perp_axis(side_face)
  const coord_a = fe_mesh_coord(a, fe, mesh)
  const is_lesser_side = side_face_is_lesser_on_perp_axis(side_face)
  coord_a == 1 && is_lesser_side || coord_a == mesh.mesh_ldims[a] && !is_lesser_side
end

import Mesh.num_boundary_sides
function num_boundary_sides(mesh::RectMesh)
  const d = mesh.space_dim
  const one = convert(Uint64, 1)
  bsides = zero(Uint64)
  for side_perp_axis=1:d
    prod = one
    for r=1:d
      prod *= (r == side_perp_axis ? 2 : mesh.mesh_ldims[r])
    end
    bsides += prod
  end
  bsides
end

import Mesh.shape_diameter_inv
shape_diameter_inv(shape::OShapeNum, mesh::RectMesh) =
  mesh.rect_diameter_inv

import Mesh.max_fe_diameter
max_fe_diameter(mesh::RectMesh) =
  mesh.rect_diameter

import Mesh.fe_interior_origin!
function fe_interior_origin!(fe::FENum, coords::Vector{R}, mesh::RectMesh)
  const d = mesh.space_dim
  for r=dim(1):d
    coords[r] = mesh.min_bounds[r] + (fe_mesh_coord(r, fe, mesh) - 1) * mesh.fe_dims[r]
  end
end


# integration functions

# Local Origins in Integration Functions
# --------------------------------------
# For each face in the mesh, the mesh must assign a local origin to be used to
# evaluate face-local functions. In this implementation, for each face F we choose
# the coordinate minimums vertex for the face, whose r^th coordinate is
#   o_r(F) = min {x_r | x in F}
# One consequence of this is that a side's local origin is the same as the interior
# origin except for the case of the greater side along a given perpendicular axis r,
# where component r of the local origin will have the constant value of that component
# for the side.
# Another consequence is that because the sides are coordinate aligned, on a side S
# perpendicular to axis r, only 0 values will be supplied in component r to functions
# which are locally defined on the side.


import Mesh.integral_face_rel_on_oshape_face
function integral_face_rel_on_oshape_face(mon::Monomial,
                                          fe_oshape::OShapeNum, face::FEFaceNum,
                                          mesh::RectMesh)
  if face == Mesh.interior_face
    Poly.integral_on_rect_at_origin(mon, mesh.fe_dims)
  else
    const a = side_face_perp_axis(face)
    dim_reduced_intgd = Poly.reduce_dim_by_fixing(a, zeroR, mon)
    Poly.integral_on_rect_at_origin(dim_reduced_intgd, mesh.fe_dims_wo_dim[a])
  end
end

import Mesh.integral_global_x_face_rel_on_fe_face
function integral_global_x_face_rel_on_fe_face(f::Function,
                                               mon::Monomial,
                                               fe::FENum, face::FEFaceNum,
                                               mesh::RectMesh)
  const d = mesh.space_dim
  const fe_x = mesh.intgd_args_work_array
  const fe_int_origin = Mesh.fe_interior_origin(fe, mesh)

  if face == Mesh.interior_face
    function ref_intgd(x::Vector{R})
      for i=1:d
        fe_x[i] = fe_int_origin[i] + x[i]
      end
      f(fe_x) * Poly.value_at(mon, x)
    end

    hcubature(ref_intgd,
              mesh.space_dim_zeros, mesh.fe_dims,
              mesh.integration_rel_err, mesh.integration_abs_err)[1]

  else # side face
    const a = side_face_perp_axis(face)
    const a_coord_of_fe_side = fe_int_origin[a] + (side_face_is_lesser_on_perp_axis(face) ? zeroR : mesh.fe_dims[a])
    const mon_dim_reduced_poly = Poly.reduce_dim_by_fixing(a, zeroR, mon)

    function ref_intgd(x::Vector{R}) # x has d-1 components
      for i=1:a-1
        fe_x[i] = fe_int_origin[i] + x[i]
      end
      fe_x[a] = a_coord_of_fe_side
      for i=a+1:d
        fe_x[i] = fe_int_origin[i] + x[i-1]
      end
      f(fe_x) * Poly.value_at(mon_dim_reduced_poly, x)
    end

    hcubature(ref_intgd,
              mesh.space_dim_less_one_zeros, mesh.fe_dims_wo_dim[a],
              mesh.integration_rel_err, mesh.integration_abs_err)[1]
  end
end

# Integrate a side-relative monomial m on its side face vs. an fe-relative vector monomial dotted with the
# outward normal.
import Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side
function integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(m::Monomial,
                                                                     q::VectorMonomial,
                                                                     fe_oshape::OShapeNum, side_face::FEFaceNum,
                                                                     mesh::RectMesh)
  const a = side_face_perp_axis(side_face)
  const qa = q[a]
  if qa == zeroR
    zeroR
  else
    const is_lesser_side = side_face_is_lesser_on_perp_axis(side_face)
    const side_fe_rel_a_coord = is_lesser_side ? zeroR : mesh.fe_dims[a]
    const qa_dim_red = Poly.reduce_dim_by_fixing(a, side_fe_rel_a_coord, qa)
    const m_dim_red = Poly.reduce_dim_by_fixing(a, zeroR, m)
    const int_m_qa = Poly.integral_on_rect_at_origin(m_dim_red * qa_dim_red, mesh.fe_dims_wo_dim[a])
    is_lesser_side ? -int_m_qa : int_m_qa
  end
end

import Mesh.integral_fe_rel_x_side_rel_on_oshape_side
function integral_fe_rel_x_side_rel_on_oshape_side(fe_mon::Monomial,
                                                   side_mon::Monomial,
                                                   fe_oshape::OShapeNum, side_face::FEFaceNum,
                                                   mesh::RectMesh)
  const a = side_face_perp_axis(side_face)
  const is_lesser_side = side_face_is_lesser_on_perp_axis(side_face)
  const side_fe_rel_a_coord = is_lesser_side ? zeroR : mesh.fe_dims[a]
  const fe_mon_dim_red = Poly.reduce_dim_by_fixing(a, side_fe_rel_a_coord, fe_mon)
  const side_mon_dim_red = Poly.reduce_dim_by_fixing(a, zeroR, side_mon)
  Poly.integral_on_rect_at_origin(fe_mon_dim_red * side_mon_dim_red, mesh.fe_dims_wo_dim[a])
end

##
##############################################



# Returns one coordinate of a finite element in the main fe/interiors mesh.
function fe_mesh_coord(r::Dim, fe::FENum, mesh::RectMesh)
  # The r^th coordinate of side n is
  #   π(r,n) = ((n − 1) mod (k_1 ··· k_r)) \ (k_1 ··· k_(r−1)) + 1
  # where k_i is the i^th component of the mesh dimensions.
  # See Rectangular_Meshes.pdf document for the derivation.
  assert(r <= mesh.space_dim, "coordinate number out of range")
  assert(fe <= mesh.num_fes, "finite element number out of range")
  mesh_coord(div(mod(fe-1, mesh.cumprods_mesh_ldims[r]), r==1 ? 1 : mesh.cumprods_mesh_ldims[r-1]) + 1)
end

function fe_mesh_coords(fe::FENum, mesh::RectMesh)
  const d = mesh.space_dim
  const coords = Array(MeshCoord, d)
  fe_mesh_coords!(fe, coords, mesh)
  coords
end

function fe_mesh_coords!(fe::FENum, coords::Array{MeshCoord,1}, mesh::RectMesh)
  const d = mesh.space_dim
  for r=dim(1):dim(d)
    coords[r] = fe_mesh_coord(r, fe, mesh)
  end
end

# Converts finite element/interior coords in the main mesh to a finite element/interior number.
function fe_with_mesh_coords(coords::Array{MeshCoord,1}, mesh::RectMesh)
  # The finite element (or interior) number for given mesh coordinates (c_1,...,c_d) is
  #   i_#(c_1,...,c_d) = 1 + sum_{i=1..d} { (c_i - 1) prod_{l=1..i-1} k_l }
  #                    = c_1 + sum_{i=2..d} { (c_i - 1) prod_{l=1..i-1} k_l }
  # where k_l is the l^th component of the mesh dimensions.
  sum = coords[1]
  for i=2:mesh.space_dim
    sum += (coords[i]-1) * mesh.cumprods_mesh_ldims[i-1]
  end
  sum
end




# side-related functions

# Find the axis which is perpendicular to the given side face.
side_face_perp_axis(side::FEFaceNum) = dim(div(side-1, 2) + 1)

# Determine whether a side face is the one with lesser axis value along its perpendicular axis.
side_face_is_lesser_on_perp_axis(side::FEFaceNum) = mod(side-1, 2) == 0

# Returns the side face with lesser axis value along the indicated axis.
lesser_side_face_perp_to_axis(a::Dim) = fefacenum(2*a - 1)

# Returns the side face with greater axis value along the indicated axis.
greater_side_face_perp_to_axis(a::Dim) = fefacenum(2*a)


# Finds the perpendicular axis for a given non-boundary side in the mesh.
function perp_axis_for_nb_side(n::NBSideNum, mesh::RectMesh)
  assert(0 < n <= mesh.num_nb_sides, "non-boundary side number out of range")
  for i=mesh.space_dim:-1:1
    if n >= mesh.first_nb_side_nums_by_perp_axis[i]
      return dim(i)
    end
  end
  error("cannot find perpendicular axis for non-boundary side number $n")
end

# Returns the geometric information for a non-boundary side in the mesh, which identifies the perpendicular
# axis for the side (and thus its orientation-specific side mesh), together with its coordinates in its
# orientation-specific mesh.
function nb_side_geom(n::NBSideNum, mesh::RectMesh)
  # The r^th coordinate of side n in the mesh of sides having the same orientation is
  #   π_s(r,n) = ((n − s_a(n)) mod Prod_{i=1..r} k_{a(n),i}) \ (Prod_{i=1..r-1} k_{a(n),i}) + 1    (r = 1,...,d)
  # where
  #   s_j is the number of the first side in the nb-side enumeration perpendicular to axis j
  #   a(n) is the axis number to which side n is perpendicular
  #   k_{j,i} is the i^th component of the dimensions of the mesh of sides perpendicular to axis j
  # See Rectangular_Meshes.pdf document for the derivation.
  const a = perp_axis_for_nb_side(n, mesh)
  const coords = Array(MeshCoord, mesh.space_dim)
  const side_mesh_rel_ix = n - mesh.first_nb_side_nums_by_perp_axis[a]
  const cumprods_side_mesh_ldims = mesh.cumprods_nb_side_mesh_ldims_by_perp_axis[a]
  # first coord is a special case because of the empty product (improper range) in the denominator
  coords[1] = mod(side_mesh_rel_ix, cumprods_side_mesh_ldims[1]) + 1
  for r=2:mesh.space_dim
    coords[r] = div(mod(side_mesh_rel_ix, cumprods_side_mesh_ldims[r]),
                    cumprods_side_mesh_ldims[r-1]) + 1
  end
  NBSideGeom(a, coords)
end

# Converts a side with a given perpendicular axis and orientation-specific side mesh coordinates
# to a non-boundary side number.
function nb_side_with_mesh_coords(coords::Array{MeshCoord,1}, perp_axis::Dim, mesh::RectMesh)
  # The enumeration number for a non-boundary side perpendicular to a given axis a, with mesh
  # coordinates (c_1,...,c_d) in its orientation-specific non-boundary side mesh, is
  #   s_{a,#}(c_1,...,c_d) = s_a + sum_{i=1..d} { (c_i - 1) prod_{l=1..i-1} k_{a,l} }
  #                        = s_a + (c_1-1) + sum_{i=2..d} { (c_i - 1) prod_{l=1..i-1} k_{a,l} }
  # where s_a is the enumeration number of the first axis-a perpendicular non-boundary side,
  # and   k_{a,l} is the l^th component of the dimensions of the mesh of axis-a perpendicular
  # non-boundary sides.
  const s_a = mesh.first_nb_side_nums_by_perp_axis[perp_axis]
  const cumprods_side_mesh_ldims = mesh.cumprods_nb_side_mesh_ldims_by_perp_axis[perp_axis]
  sum = s_a + coords[1] - 1
  for i=2:mesh.space_dim
    sum += (coords[i]-1) * cumprods_side_mesh_ldims[i-1]
  end
  sum
end

# Exporting

function exportAsGmshSurface(ios::IO,  mesh::RectMesh)
  if mesh.space_dim != 2 
    error("Gmsh output for rectangle meshes is currently only supported for the 2d case.")
  end

  typealias PointNum Uint64
  typealias LineNum Int64

  const pointnums_by_mcoords = sizehint(Dict{(MeshCoord, MeshCoord), PointNum}(), uint64(mesh.num_fes))
  const linenums_by_endptnums = sizehint(Dict{(PointNum, PointNum), LineNum}(), uint64(2*mesh.num_fes))

  const incrementor() = let next = uint64(1); () -> (next += 1) - 1 end
  const next_pointnum = incrementor()
  const next_linenum  = incrementor()

  function register_point(pt_mcoords::(MeshCoord, MeshCoord))
    const existing_pointnum = get(pointnums_by_mcoords, pt_mcoords, uint64(0))
    if existing_pointnum == 0
      # Register the new point number with its mesh coordinates.
      const pointnum = next_pointnum()
      pointnums_by_mcoords[pt_mcoords] = pointnum
      # Write the Gmsh point declaration.
      const pt_1 = mesh.min_bounds[1] + (pt_mcoords[1] - 1) * mesh.fe_dims[1]
      const pt_2 = mesh.min_bounds[2] + (pt_mcoords[2] - 1) * mesh.fe_dims[2]
      @printf(ios, "Point(%u) = {%.15le, %.15le, 0.0, 1.0};\n", pointnum, pt_1, pt_2)
      pointnum
    else
      existing_pointnum
    end
  end

  # Fetch the negative of a registered line number for a pair of points if the points have already been registered in
  # the reverse order (which should be the only case in which the line is found at all), else register and return a
  # (positive) new number for the line.
  function get_linenum(ptnum_1::PointNum, ptnum_2::PointNum)
    const minus_existing_linenum = -get(linenums_by_endptnums, (ptnum_2, ptnum_1), uint64(0))
    if minus_existing_linenum != 0
      minus_existing_linenum
    else
      # Register new line.
      const linenum = next_linenum()
      linenums_by_endptnums[(ptnum_1, ptnum_2)] = linenum
      linenum
    end
  end

  # pre-allocated work arrays for the loop
  const fe_mcoords = Array(MeshCoord, mesh.space_dim)

  for fe=Mesh.fenum(1):mesh.num_fes
    fe_mesh_coords!(fe, fe_mcoords, mesh)

    # lower left vertex
    const ll_pointnum = register_point((fe_mcoords[1], fe_mcoords[2]))

    # lower right vertex
    const lr_pointnum = register_point((fe_mcoords[1]+1, fe_mcoords[2]))

    # lower left to lower right line
    const ll_lr_linenum = get_linenum(ll_pointnum, lr_pointnum)
    if ll_lr_linenum > 0 # new line: declare in Gmsh
      @printf(ios, "Line(%d) = {%u, %u};\n", ll_lr_linenum, ll_pointnum, lr_pointnum)
    end

    # upper right vertex
    const ur_pointnum = register_point((fe_mcoords[1]+1, fe_mcoords[2]+1))

    # lower right to upper right line
    const lr_ur_linenum = get_linenum(lr_pointnum, ur_pointnum)
    if lr_ur_linenum > 0 # new line: declare in Gmsh
      @printf(ios, "Line(%d) = {%u, %u};\n", lr_ur_linenum, lr_pointnum, ur_pointnum)
    end

    # upper left vertex
    const ul_pointnum = register_point((fe_mcoords[1], fe_mcoords[2]+1))

    # upper right to upper left line
    const ur_ul_linenum = get_linenum(ur_pointnum, ul_pointnum)
    if ur_ul_linenum > 0 # new line: declare in Gmsh
      @printf(ios, "Line(%d) = {%u, %u};\n", ur_ul_linenum, ur_pointnum, ul_pointnum)
    end

    # upper left to lower left line
    const ul_ll_linenum = get_linenum(ul_pointnum, ll_pointnum)
    if ul_ll_linenum > 0 # new line: declare in Gmsh
      @printf(ios, "Line(%d) = {%u, %u};\n", ul_ll_linenum, ul_pointnum, ll_pointnum)
    end

    # Write Gmsh line loop representing this element.
    @printf(ios, "Line Loop(%u) = {%d, %d, %d, %d};\n", fe, ll_lr_linenum, lr_ur_linenum, ur_ul_linenum, ul_ll_linenum)

    # Write Gmsh surface representing this element constructed from the above line loop.
    @printf(ios, "Plane Surface(%u) = {%u};\n", fe, fe)

    # Tag the physical region with this finite element number.
    @printf(ios, "Physical Line(%u) = {%d, %d, %d, %d};\n\n\n", fe, ll_lr_linenum, lr_ur_linenum, ur_ul_linenum, ul_ll_linenum)
  end

  flush(ios)
end


end # end of module
