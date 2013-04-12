module RMesh3
export RectMesh3,
       MeshCoord, mesh_coord,
       NBSideInfo, nb_side_info,
       FEInfo, fe_info

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FERelFace, Mesh.fe_face, Mesh.AbstractMesh, Mesh.NBSideInclusions
import Poly, Poly.Monomial, Poly.VectorMonomial, Poly.Polynomial
import Cubature.hcubature

# Abbreviations
# nb: non-boundary
# A_nb_side: non-boundary side perpendicular to the A-axis (A=x,y,z)
#

# type for mesh coords (columns, rows, stacks)
typealias MeshCoord Uint64
mesh_coord(i::Integer) = convert(MeshCoord, i)

const default_integration_rel_err = 10e-8
const default_integration_abs_err = 10e-8


type RectMesh3 <: AbstractMesh

  min_x::R
  min_y::R
  min_z::R

  max_x::R
  max_y::R
  max_z::R

  cols::MeshCoord
  rows::MeshCoord
  stacks::MeshCoord

  # Computed values

  fe_dims::Array{R,1}
  fe_dims_wo_1::Array{R,1}
  fe_dims_wo_2::Array{R,1}
  fe_dims_wo_3::Array{R,1}

  num_elements::FENum

  num_x_nb_sides::NBSideNum
  num_y_nb_sides::NBSideNum
  num_z_nb_sides::NBSideNum
  num_nb_sides::NBSideNum

  sidenum_first_x_nb_side::NBSideNum
  sidenum_first_y_nb_side::NBSideNum
  sidenum_first_z_nb_side::NBSideNum

  num_fes_per_stack::NBSideNum
  num_x_nb_sides_per_stack::NBSideNum
  num_y_nb_sides_per_stack::NBSideNum
  num_z_nb_sides_per_stack::NBSideNum

  # integration support members
  intgd_args_work_array::Vector{R}
  ref_fe_min_bounds::Vector{R}
  ref_fe_min_bounds_short::Vector{R}
  integration_rel_err::R
  integration_abs_err::R


  function RectMesh3(min_coords::(R,R,R), max_coords::(R,R,R),
                     ncols::Integer, nrows::Integer, nstacks::Integer,
                     integration_rel_err::R, integration_abs_err::R)
    assert(ncols > 0 && nrows > 0 && nstacks > 0, "positive columns, rows and stacks required")
    assert(max_coords[1] > min_coords[1] && max_coords[2] > min_coords[2] && max_coords[3] > min_coords[3], "min exceeds max")

    const num_x_nb_sides = (ncols-1) * nrows * nstacks
    const num_y_nb_sides = ncols * (nrows-1) * nstacks
    const num_z_nb_sides = ncols * nrows * (nstacks-1)

    const sidenum_first_x_nb_side = 1
    const sidenum_first_y_nb_side = sidenum_first_x_nb_side + num_x_nb_sides
    const sidenum_first_z_nb_side = sidenum_first_y_nb_side + num_y_nb_sides

    const fes_per_stack = ncols * nrows

    const fe_dims = [(max_coords[1] - min_coords[1]) / ncols,
                     (max_coords[2] - min_coords[2]) / nrows,
                     (max_coords[3] - min_coords[3]) / nstacks]
    const fe_dims_wo_1 = fe_dims[2:3]
    const fe_dims_wo_2 = fe_dims[[1,3]]
    const fe_dims_wo_3 = fe_dims[1:2]

    new(min_coords[1], min_coords[2], min_coords[3],
        max_coords[1], max_coords[2], max_coords[3],
        convert(MeshCoord, ncols), convert(MeshCoord, nrows), convert(MeshCoord, nstacks),
        fe_dims, fe_dims_wo_1, fe_dims_wo_2, fe_dims_wo_3,
        nrows * ncols * nstacks, # total fes
        num_x_nb_sides,
        num_y_nb_sides,
        num_z_nb_sides,
        num_x_nb_sides + num_y_nb_sides + num_z_nb_sides, # num_nb_sides
        sidenum_first_x_nb_side,
        sidenum_first_y_nb_side,
        sidenum_first_z_nb_side,
        fes_per_stack,
        (ncols-1) * nrows, # num_x_nb_sides_per_stack
        ncols * (nrows-1), # num_y_nb_sides_per_stack
        fes_per_stack,     # num_z_nb_sides_per_stack
        Array(R, 3), # integrand args work array
        zeros(R, 3), # ref fe min bounds
        zeros(R, 2), # ref fe min bounds, short
        integration_rel_err,
        integration_abs_err)
  end # RectMesh3 constructor
end # type RectMesh3

RectMesh3(min_coords::(R,R,R), max_coords::(R,R,R),
          ncols::Integer, nrows::Integer, nstacks::Integer) =
  RectMesh3(min_coords, max_coords, ncols, nrows, nstacks, default_integration_rel_err, default_integration_abs_err)


#############################################
# Implement functions required of all meshes.

import Mesh.space_dim
space_dim(m::RectMesh3) = dim(3)

const one_monomial = Monomial(0,0,0)
import Mesh.one_mon
one_mon(m::RectMesh3) = one_monomial

import Mesh.num_fes
num_fes(mesh::RectMesh3) = mesh.num_elements

import Mesh.num_nb_sides
num_nb_sides(mesh::RectMesh3) = mesh.num_nb_sides

import Mesh.num_side_faces_per_fe
num_side_faces_per_fe(mesh::RectMesh3) = 6

import Mesh.dependent_dim_for_nb_side
dependent_dim_for_nb_side(i::NBSideNum, mesh::RectMesh3) =
  is_x_nb_side(i, mesh) ? 1 : is_y_nb_side(i, mesh) ? 2 : 3

import Mesh.dependent_dim_for_ref_side_face
dependent_dim_for_ref_side_face(side_face::FERelFace, mesh::RectMesh3) = side_face_perp_axis(side_face)

import Mesh.fe_inclusions_of_nb_side!
function fe_inclusions_of_nb_side!(i::NBSideNum, mesh::RectMesh3, fe_incls::NBSideInclusions)
  fe_incls.nb_side_num = i
  if is_x_nb_side(i, mesh)
    const side_type_rel_ix = i - mesh.sidenum_first_x_nb_side
    const sides_per_stack = mesh.num_x_nb_sides_per_stack
    const stack_ix =     div(side_type_rel_ix, sides_per_stack)
    const stack_rel_ix = mod(side_type_rel_ix, sides_per_stack)
    const side_cols = mesh.cols-1
    const side_row_ix = div(stack_rel_ix, side_cols)
    const side_col_ix = mod(stack_rel_ix, side_cols)
    # now find the fe's which adjoin at the side
    const fe1 = stack_ix * mesh.num_fes_per_stack + side_row_ix * mesh.cols + side_col_ix + 1
    fe_incls.fe1 = fe1
    fe_incls.face_in_fe1 = x_max_face # not backwards, max face in the lesser fe in this direction
    fe_incls.fe2 = fe1 + 1
    fe_incls.face_in_fe2 = x_min_face # min face of greater fe
  elseif is_y_nb_side(i, mesh)
    const side_type_rel_ix = i - mesh.sidenum_first_y_nb_side
    const sides_per_stack = mesh.num_y_nb_sides_per_stack
    const stack_ix =     div(side_type_rel_ix, sides_per_stack)
    const stack_rel_ix = mod(side_type_rel_ix, sides_per_stack)
    const cols = mesh.cols # same number of cols for sides and fes in this case
    const side_row_ix = div(stack_rel_ix, cols)
    const side_col_ix = mod(stack_rel_ix, cols)
    # now find the fe's which adjoin at the side
    const fe1 = stack_ix * mesh.num_fes_per_stack + side_row_ix * cols + side_col_ix + 1
    fe_incls.fe1 = fe1
    fe_incls.face_in_fe1 = y_max_face # not backwards, max face in the lesser fe in this direction
    fe_incls.fe2 = fe1 + cols
    fe_incls.face_in_fe2 = y_min_face # min face of greater fe
  elseif is_z_nb_side(i, mesh)
    const side_type_rel_ix = i - mesh.sidenum_first_z_nb_side
    const sides_per_stack = mesh.num_z_nb_sides_per_stack
    const stack_ix =     div(side_type_rel_ix, sides_per_stack)
    const stack_rel_ix = mod(side_type_rel_ix, sides_per_stack)
    const cols = mesh.cols # same number of cols for sides and fes in this case
    const side_row_ix = div(stack_rel_ix, cols)
    const side_col_ix = mod(stack_rel_ix, cols)
    # now find the fe's which adjoin at the side
    const fe1 = stack_ix * mesh.num_fes_per_stack + side_row_ix * cols + side_col_ix + 1
    fe_incls.fe1 = fe1
    fe_incls.face_in_fe1 = z_max_face # not backwards, max face in the lesser fe in this direction
    fe_incls.fe2 = fe1 + mesh.num_fes_per_stack
    fe_incls.face_in_fe2 = z_min_face # min face of greater fe
  else
    error("invalid side number")
  end
end

import Mesh.is_boundary_side
function is_boundary_side(fe::FENum, face::FERelFace, mesh::RectMesh3)
  if face == x_min_face
    fe_col(fe, mesh) == 1
  elseif face == x_max_face
    fe_col(fe, mesh) == mesh.cols
  elseif face == y_min_face
    fe_row(fe, mesh) == 1
  elseif face == y_max_face
    fe_row(fe, mesh) == mesh.rows
  elseif face == z_min_face
    fe_stack(fe, mesh) == 1
  elseif face == z_max_face
    fe_stack(fe, mesh) == mesh.stacks
  else
    error("invalid face: $face")
  end
end

import Mesh.integral_face_rel_on_face
integral_face_rel_on_face(mon::Monomial, face::FERelFace, mesh::RectMesh3) =
  if face == Mesh.interior_face
    Poly.integral_on_rect_at_origin(mon, mesh.fe_dims)
  elseif face == x_min_face
    if mon.exps[1] != 0
      zeroR
    else
      dim_reduced = Poly.reduce_dim_by_fixing(dim(1), zeroR, mon)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_1)
    end
  elseif face == x_max_face
    dim_reduced = Poly.reduce_dim_by_fixing(dim(1), mesh.fe_dims[1], mon)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_1)
  elseif face == y_min_face
    if mon.exps[2] != 0
      zeroR
    else
      dim_reduced = Poly.reduce_dim_by_fixing(dim(2), zeroR, mon)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_2)
    end
  elseif face == y_max_face
    dim_reduced = Poly.reduce_dim_by_fixing(dim(2), mesh.fe_dims[2], mon)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_2)
  elseif face == z_min_face
    if mon.exps[3] != 0
      zeroR
    else
      dim_reduced = Poly.reduce_dim_by_fixing(dim(3), zeroR, mon)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_3)
    end
  elseif face == z_max_face
    dim_reduced = Poly.reduce_dim_by_fixing(dim(3), mesh.fe_dims[3], mon)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_dims_wo_3)
  else
    error("invalid face: $face")
  end



import Mesh.integral_on_ref_fe_side_vs_outward_normal
function integral_on_ref_fe_side_vs_outward_normal(vm::VectorMonomial, face::FERelFace, mesh::RectMesh3)
  if face == x_min_face
    -integral_face_rel_on_face(vm[1], face, mesh)
  elseif face == x_max_face
    integral_face_rel_on_face(vm[1], face, mesh)
  elseif face == y_min_face
    -integral_face_rel_on_face(vm[2], face, mesh)
  elseif face == y_max_face
    integral_face_rel_on_face(vm[2], face, mesh)
  elseif face == z_min_face
    -integral_face_rel_on_face(vm[3], face, mesh)
  elseif face == z_max_face
     integral_face_rel_on_face(vm[3], face, mesh)
  else
    error("invalid face: $face")
  end
end

# TODO: unit tests
import Mesh.integral_global_x_face_rel_on_fe_face
function integral_global_x_face_rel_on_fe_face(f::Function, mon::Monomial, fe::FENum, face::FERelFace, mesh::RectMesh3)
  const d = 3
  const fe_local_origin = fe_coords(fe, mesh)
  const fe_x = mesh.intgd_args_work_array
  if face == Mesh.interior_face
    function ref_intgd(x::Vector{R})
      for r=1:d
        fe_x[r] = fe_local_origin[r] + x[r]
      end
      f(fe_x) * Poly.monomial_value(mon, x)
    end
    hcubature(ref_intgd, mesh.ref_fe_min_bounds, mesh.fe_dims, mesh.integration_rel_err, mesh.integration_abs_err)[1]
  else # side face
    # perp axis for the side
    const a = face == x_min_face || face == x_max_face ? 1 : face == y_min_face || face == y_max_face ? 2 : 3
    const a_coord_of_ref_side = face == x_min_face || face == y_min_face || face == z_min_face ? zeroR : mesh.fe_dims[a]
    const a_coord_of_fe_side = fe_local_origin[a] + a_coord_of_ref_side
    const dim_reduced_poly = Poly.reduce_dim_by_fixing(dim(a), a_coord_of_ref_side, mon)
    function ref_intgd(x::Vector{R}) # x has 1 component
      for i=1:a-1
        fe_x[i] = fe_local_origin[i] + x[i]
      end
      fe_x[a] = a_coord_of_fe_side
      for i=a+1:d
        fe_x[i] = fe_local_origin[i] + x[i-1]
      end
      f(fe_x) * Poly.polynomial_value(dim_reduced_poly, x)
    end
    const fe_dims_wo_a = a == 1 ? mesh.fe_dims_wo_1 : a == 2 ? mesh.fe_dims_wo_2 : mesh.fe_dims_wo_3
    hcubature(ref_intgd, mesh.ref_fe_min_bounds_short, fe_dims_wo_a, mesh.integration_rel_err, mesh.integration_abs_err)[1]
  end
end


#
##############################################


# functions specific to rectangular meshes

is_x_nb_side(i::NBSideNum, mesh::RectMesh3) = mesh.sidenum_first_x_nb_side <= i < mesh.sidenum_first_y_nb_side
is_y_nb_side(i::NBSideNum, mesh::RectMesh3) = mesh.sidenum_first_y_nb_side <= i < mesh.sidenum_first_z_nb_side
is_z_nb_side(i::NBSideNum, mesh::RectMesh3) = mesh.sidenum_first_z_nb_side <= i <= mesh.num_nb_sides

side_face_perp_axis(side_face::FERelFace) =
  face == x_min_face || face == x_max_face ? 1 : face == y_min_face || face == y_max_face ? 2 : 3

# face numbers
# Only side faces are defined here, the interior_face is defined in the Mesh module.
const x_min_face = fe_face(1)
const x_max_face = fe_face(2)
const y_min_face = fe_face(3)
const y_max_face = fe_face(4)
const z_min_face = fe_face(5)
const z_max_face = fe_face(6)
const num_sides = 0x06

function fe_row_ix(fe_ix::FENum, mesh::RectMesh3)
  const stack_rel_ix = mod(fe_ix, mesh.num_fes_per_stack)
  div(stack_rel_ix, mesh.cols)
end

function fe_col_ix(fe_ix::FENum, mesh::RectMesh3)
  const stack_rel_ix = mod(fe_ix, mesh.num_fes_per_stack)
  mod(stack_rel_ix, mesh.cols)
end

function fe_stack_ix(fe_ix::FENum, mesh::RectMesh3)
  div(fe_ix, mesh.num_fes_per_stack)
end

fe_row(fe::FENum, mesh::RectMesh3) = fe_row_ix(fe-1, mesh) + 1
fe_col(fe::FENum, mesh::RectMesh3) = fe_col_ix(fe-1, mesh) + 1
fe_stack(fe::FENum, mesh::RectMesh3) = fe_stack_ix(fe-1, mesh) + 1


# Fill the passed 3 element array with the coordinates of the finite element corner with range-minimum coordinates.
function fe_coords!(fe::FENum, mesh::RectMesh3, coords::Vector{R})
  const fe_ix = fe - 1
  coords[1] = mesh.min_x + fe_col_ix(fe_ix, mesh)   * mesh.fe_dims[1]
  coords[2] = mesh.min_y + fe_row_ix(fe_ix, mesh)   * mesh.fe_dims[2]
  coords[3] = mesh.min_z + fe_stack_ix(fe_ix, mesh) * mesh.fe_dims[3]
end

# Functional variant of the above.
fe_coords(fe::FENum, mesh::RectMesh3) =
  let a = zeros(R,3)
    fe_coords!(fe, mesh, a)
    a
  end


###############################################
# side and fe information retrieval for testing

type FEInfo
  fe_num::FENum
  col::MeshCoord
  row::MeshCoord
  stack::MeshCoord
  coords::Array{R,1}
  function FEInfo(fe::FENum, col::Integer, row::Integer, stack::Integer, coords::Array{R,1})
    new(Mesh.fe_num(fe),
        convert(MeshCoord,col), convert(MeshCoord,row), convert(MeshCoord,stack),
        coords)
  end
end

import Base.isequal
isequal(fei1::FEInfo, fei2::FEInfo) =
  fei1.fe_num == fei2.fe_num &&
  fei1.col    == fei2.col &&
  fei1.row    == fei2.row &&
  fei1.stack  == fei2.stack &&
  isequal(fei1.coords, fei2.coords)

import Base.hash
hash(fei::FEInfo) = fei.fe_num + 3*fei.col + 5*fe.row + 7*fei.stack + 11*hash(fei.coords)

fe_info(fe::FENum, mesh::RectMesh3) =
  FEInfo(fe,
         fe_col(fe, mesh),
         fe_row(fe, mesh),
         fe_stack(fe, mesh),
         fe_coords(fe, mesh))

import Base.string, Base.show, Base.print
function string(fei::FEInfo)
  "FEInfo: fe: $(dec(fei.fe_num)), col:$(dec(fei.col)), row:$(dec(fei.row)), stack:$(dec(fei.stack)), coords:$(fei.coords)"
end
print(io::IO, fei::FEInfo) = print(io, string(fei))
show(io::IO, fei::FEInfo) = print(io, fei)


type NBSideInfo
  side_num::NBSideNum
  perp_axis::Char
  lesser_adjoining_fe::FEInfo
  greater_adjoining_fe::FEInfo
end

function nb_side_info(side_num::NBSideNum, mesh::RectMesh3)
  incls = Mesh.fe_inclusions_of_nb_side(side_num, mesh)
  perp_to_axis = is_x_nb_side(side_num, mesh) ? 'x' : is_y_nb_side(side_num, mesh) ? 'y' : 'z'

  NBSideInfo(side_num, perp_to_axis, fe_info(incls.fe1, mesh), fe_info(incls.fe2, mesh))
end

isequal(si1::NBSideInfo, si2::NBSideInfo) =
  si1.side_num  == si2.side_num &&
  si1.perp_axis == si2.perp_axis &&
  isequal(si1.lesser_adjoining_fe, si2.lesser_adjoining_fe) &&
  isequal(si1.greater_adjoining_fe, si2.greater_adjoining_fe)

import Base.hash
hash(si::NBSideInfo) = si.side_num + 3*si.perp_axis + 5*hash(si.lesser_adjoining_fe) + 7*hash(si.greater_adjoining_fe)

# string representation
function string(si::NBSideInfo)
  "NBSideInfo:  side:$(dec(si.side_num))\n  perp to axis: $(si.perp_axis)\n  lesser_fe:  $(si.lesser_adjoining_fe)\n  greater_fe: $(si.greater_adjoining_fe)"
end
print(io::IO, si::NBSideInfo) = print(io, string(si))
show(io::IO, si::NBSideInfo) = print(io, si)


end # end of module
