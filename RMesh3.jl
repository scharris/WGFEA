module RMesh3
export RectMesh3,
       SideInfo, side_info,
       FEInfo, fe_info

using Common
import Mesh, Mesh.FENum, Mesh.SideNum, Mesh.FEFace, Mesh.fe_face, Mesh.AbstractMesh, Mesh.NBSideInclusions
import Poly, Poly.Monomial, Poly.VectorMonomial

# Abbreviations
# nb: non-boundary
# A_nb_side: non-boundary side perpendicular to the A-axis (A=x,y,z)
#

# type for mesh coords (columns, rows, stacks)
typealias MeshAxisVal Uint64

type RectMesh3 <: AbstractMesh

  min_x::R
  min_y::R
  min_z::R

  max_x::R
  max_y::R
  max_z::R

  cols::MeshAxisVal
  rows::MeshAxisVal
  stacks::MeshAxisVal

  # Computed values

  fe_width::R
  fe_height::R
  fe_depth::R

  num_elements::FENum

  num_x_nb_sides::SideNum
  num_y_nb_sides::SideNum
  num_z_nb_sides::SideNum
  num_nb_sides::SideNum

  sidenum_first_x_nb_side::SideNum
  sidenum_first_y_nb_side::SideNum
  sidenum_first_z_nb_side::SideNum

  num_fes_per_stack::SideNum
  num_x_nb_sides_per_stack::SideNum
  num_y_nb_sides_per_stack::SideNum
  num_z_nb_sides_per_stack::SideNum

  function RectMesh3(min_coords::(R,R,R), max_coords::(R,R,R), ncols::Integer, nrows::Integer, nstacks::Integer)
    assert(ncols > 0 && nrows > 0 && nstacks > 0, "positive columns, rows and stacks required")
    assert(max_coords[1] > min_coords[1] && max_coords[2] > min_coords[2] && max_coords[3] > min_coords[3], "min exceeds max")

    const num_x_nb_sides = (ncols-1) * nrows * nstacks
    const num_y_nb_sides = ncols * (nrows-1) * nstacks
    const num_z_nb_sides = ncols * nrows * (nstacks-1)

    const sidenum_first_x_nb_side = 1
    const sidenum_first_y_nb_side = sidenum_first_x_nb_side + num_x_nb_sides
    const sidenum_first_z_nb_side = sidenum_first_y_nb_side + num_y_nb_sides

    const fes_per_stack = ncols * nrows

    new(min_coords[1], min_coords[2], min_coords[3],
        max_coords[1], max_coords[2], max_coords[3],
        convert(MeshAxisVal, ncols), convert(MeshAxisVal, nrows), convert(MeshAxisVal, nstacks),
        (max_coords[1] - min_coords[1]) / ncols,   # fe width
        (max_coords[2] - min_coords[2]) / nrows,   # fe height
        (max_coords[3] - min_coords[3]) / nstacks, # fe depth
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
        fes_per_stack)     # num_z_nb_sides_per_stack
  end # RectMesh3 constructor
end # type RectMesh3

## The following three functions compute 0-based row, col, and stack indexes in
## the grid consisting of items of the same type only (e.g. faces of a
## particular orientation), for a given item index relative to the first of the
## same type, when the items are enumerated by column first (changing fastest),
## then more significantly by row, and most significantly by stack. The number
## of columns and rows are passed in because they depend on the type of item
## being enumerated, e.g. nonboundary sides of constant x only have mesh.cols-1
## columns per row.
#function mesh_item_row_ix(item_ix::Integer, num_item_cols::Integer, num_item_rows::Integer)
#  const items_per_stack = num_item_rows * num_item_cols
#  const stack_rel_ix = mod(item_ix, items_per_stack)
#  div(stack_rel_ix, num_item_cols)
#end
#function mesh_item_col_ix(item_ix::Integer, num_item_cols::Integer, num_item_rows::Integer)
#  const items_per_stack = num_item_rows * num_item_cols
#  const stack_rel_ix = mod(item_ix, items_per_stack)
#  mod(stack_rel_ix, num_item_cols)
#end
#function mesh_item_stack_ix(item_ix::Integer, num_item_cols::Integer, num_item_rows::Integer)
#  const items_per_stack = num_item_rows * num_item_cols
#  div(item_ix, items_per_stack)
#end


# ------------------------------------------
# Implement functions required of all meshes.

import Mesh.space_dim
space_dim(m::RectMesh3) = dim(3)

import Mesh.num_fes
num_fes(mesh::RectMesh3) = mesh.num_elements

import Mesh.num_nb_sides
num_nb_sides(mesh::RectMesh3) = mesh.num_nb_sides

import Mesh.num_side_faces_per_fe
num_side_faces_per_fe(mesh::RectMesh3) = 6

import Mesh.fe_inclusions_of_nb_side!
function fe_inclusions_of_nb_side!(i::SideNum, mesh::RectMesh3, fe_incls::NBSideInclusions)
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

import Mesh.integral_on_ref_fe_interior
integral_on_ref_fe_interior(mon::Monomial, mesh::RectMesh3) =
  Poly.integral_on_rect_at_origin(mon, mesh.fe_width, mesh.fe_height, mesh.fe_depth)

import Mesh.integral_on_ref_fe_side_vs_outward_normal
function integral_on_ref_fe_side_vs_outward_normal(vm::VectorMonomial, face::FEFace, mesh::RectMesh3)
  if face == x_min_face
    # On the x min face, vm . n = -vm[1], and dim 1 is 0.
    # Recognize trivially zero cases early to avoid unnecessary allocations.
    if vm.mon_pos != 1 || vm[1].exps[1] != 0
      zeroR
    else
      const vm_dot_n = -vm[1]
      # Fix dim 1 and integrate the function of the remaining dimensions.
      dim_reduced = Poly.reduce_dim_by_fixing(dim(1), zeroR, vm_dot_n)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_height, mesh.depth)
    end
  elseif face == x_max_face
    # On the x max face, vm . n = vm[1], and dim 1 is constantly fe_width.
    const vm_dot_n = vm[1]
    # Fix dim 1 and integrate over the remaining non-fixed dimensions.
    dim_reduced = Poly.reduce_dim_by_fixing(dim(1), mesh.fe_width, vm_dot_n)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_height, mesh.depth)
  elseif face == y_min_face
    # On the y min face, vm . n = -vm[2], and dim 2 is constantly 0.
    if vm.mon_pos != 2 || vm[2].exps[2] != 0
      zeroR
    else
      const vm_dot_n = -vm[2]
      # Fix dim 2 and integrate over the remaining non-fixed dimensions.
      dim_reduced = Poly.reduce_dim_by_fixing(dim(2), zeroR, vm_dot_n)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_width, mesh.fe_depth)
    end
  elseif face == y_max_face
    # On the y max face, vm . n = vm[2], and dim 2 is constantly fe_height.
    const vm_dot_n = vm[2]
    # Fix dim 2 and integrate over the remaining non-fixed dimensions.
    dim_reduced = Poly.reduce_dim_by_fixing(dim(2), mesh.fe_height, vm_dot_n)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_width, mesh.fe_depth)
  elseif face == z_min_face
    # On the z min face, vm . n = -vm[3], and dim 3 is constantly 0.
    if vm.mon_pos != 3 || vm[3].exps[3] != 0
      zeroR
    else
      const vm_dot_n = -vm[3]
      # Fix dim 3 and integrate over the remaining non-fixed dimensions.
      dim_reduced = Poly.reduce_dim_by_fixing(dim(3), zeroR, vm_dot_n)
      Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_width, mesh.fe_height)
    end
  elseif face == z_max_face
    # On the z max face, vm . n = vm[3], and dim 3 is constantly fe_depth.
    const vm_dot_n = vm[3]
    # Fix dim 3 and integrate over the remaining non-fixed dimensions.
    dim_reduced = Poly.reduce_dim_by_fixing(dim(3), mesh.fe_depth, vm_dot_n)
    Poly.integral_on_rect_at_origin(dim_reduced, mesh.fe_width, mesh.fe_height)
  else
    error("invalid face: $face")
  end
end

#
# ------------------------------------------

# functions specific to rectangular meshes

is_x_nb_side(i::SideNum, mesh::RectMesh3) = mesh.sidenum_first_x_nb_side <= i < mesh.sidenum_first_y_nb_side
is_y_nb_side(i::SideNum, mesh::RectMesh3) = mesh.sidenum_first_y_nb_side <= i < mesh.sidenum_first_z_nb_side
is_z_nb_side(i::SideNum, mesh::RectMesh3) = mesh.sidenum_first_z_nb_side <= i <= mesh.num_nb_sides

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


# Fill the passed 3 element array with the coordinates of corner with range-minimum coordinates.
function fe_coords!(fe::FENum, mesh::RectMesh3, coords::Vector{R})
  const fe_ix = fe - 1
  coords[1] = mesh.min_x + fe_col_ix(fe_ix, mesh)   * mesh.fe_width
  coords[2] = mesh.min_y + fe_row_ix(fe_ix, mesh)   * mesh.fe_height
  coords[3] = mesh.min_z + fe_stack_ix(fe_ix, mesh) * mesh.fe_depth
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
  col::MeshAxisVal
  row::MeshAxisVal
  stack::MeshAxisVal
  coords::Array{R,1}
  function FEInfo(fe::FENum, col::Integer, row::Integer, stack::Integer, coords::Array{R,1})
    new(Mesh.fe_num(fe),
        convert(MeshAxisVal,col), convert(MeshAxisVal,row), convert(MeshAxisVal,stack),
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


type SideInfo
  side_num::SideNum
  perp_axis::Char
  lesser_adjoining_fe::FEInfo
  greater_adjoining_fe::FEInfo
end

function side_info(side_num::SideNum, mesh::RectMesh3)
  incls = Mesh.fe_inclusions_of_nb_side(side_num, mesh)
  perp_to_axis = is_x_nb_side(side_num, mesh) ? 'x' : is_y_nb_side(side_num, mesh) ? 'y' : 'z'

  SideInfo(side_num, perp_to_axis, fe_info(incls.fe1, mesh), fe_info(incls.fe2, mesh))
end

isequal(si1::SideInfo, si2::SideInfo) =
  si1.side_num  == si2.side_num &&
  si1.perp_axis == si2.perp_axis &&
  isequal(si1.lesser_adjoining_fe, si2.lesser_adjoining_fe) &&
  isequal(si1.greater_adjoining_fe, si2.greater_adjoining_fe)

import Base.hash
hash(si::SideInfo) = si.side_num + 3*si.perp_axis + 5*hash(si.lesser_adjoining_fe) + 7*hash(si.greater_adjoining_fe)

# string representation
function string(si::SideInfo)
  "SideInfo:  side:$(dec(si.side_num))\n  perp to axis: $(si.perp_axis)\n  lesser_fe:  $(si.lesser_adjoining_fe)\n  greater_fe: $(si.greater_adjoining_fe)"
end
print(io::IO, si::SideInfo) = print(io, string(si))
show(io::IO, si::SideInfo) = print(io, si)


end # end of module
