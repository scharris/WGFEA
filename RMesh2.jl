module RMesh2
export RectMesh2

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFace, Mesh.fe_face, Mesh.AbstractMesh, Mesh.NBSideInclusions
import Poly, Poly.Monomial, Poly.VectorMonomial, Poly.Polynomial
import Cubature.hcubature, Cubature.hquadrature

# Ordering of Faces
# -----------------
# Sides and interiors are numbered separately. All sides are numbered together, with
# vertical non-boundary sides first, followed by the horizontal non-boundary sides.
# ==interiors numbering==     =================== sides numbering ========================
# |<    1 ... rc      >|      |<      1 ...  r(c-1)         r(c-1)+1 ...  r(c-1)+(r-1)c >|
#
# With r and c being the number of rows and columns of finite elements in the mesh,
# the finite element faces are ordered as shown below:
#
# interiors ordering          vertical nb sides ordering    horizontal nb sides ordering
# ----------------------      --------------------------    ----------------------------
# |...            | rc |      |...          | r(c-1)|  |    |     |   ...   |          |
#         ...                            ...                ----------------------------
# ----------------------      --------------------------    |...            | (r-1)c^  |
# | c+1 | c+2 |...| 2c |      |  c|  c+1|...| 2(c-1)|  |              ...
# ----------------------      --------------------------    ----------------------------
# |  1  |  2  |...|  c |      |  1|    2|...|  (c-1)|  |    | 1^  |   ...   |    c^    |
# ----------------------      --------------------------    ----------------------------
# (^also fe numbering)                                      (^add # vert nb sides to get side #'s)


# Reference Finite Element
# ----------------------------------------------------
# Operations involving polynomials are performed on a reference finite element with
# origin at the lower left corner of the element.  Such results then apply to a general
# finite element for the same polynomials pre-composed with the translation function
# (x,y) -> (x-x0, y-y0), where (x0,y0) is the lower left corner of the element.


typealias MeshCoord Uint64

const default_integration_rel_err = 10e-8
const default_integration_abs_err = 10e-8

type RectMesh2 <: AbstractMesh
  bottom_left_x::R
  bottom_left_y::R
  top_right_x::R
  top_right_y::R
  rows::MeshCoord
  cols::MeshCoord
  # Computed values
  num_elements::FENum
  num_nb_sides::NBSideNum
  num_nb_vert_sides::NBSideNum
  num_nb_horz_sides::NBSideNum
  sidenum_first_nb_horz_side::NBSideNum

  fe_dims::Array{R,1}
  fe_dims_wo_1::Array{R,1}
  fe_dims_wo_2::Array{R,1}

  # integration support members
  intgd_args_work_array::Vector{R}
  ref_fe_min_bounds::Vector{R}
  ref_fe_min_bounds_short::Vector{R}
  integration_rel_err::R
  integration_abs_err::R

  function RectMesh2(bottom_left::(R,R), top_right::(R,R),
                     nrows::Integer, ncols::Integer,
                     integration_rel_err::R,
                     integration_abs_err::R)
    assert(nrows > 0 && ncols > 0, "positive rows and columns required")
    assert(top_right[1] > bottom_left[1] && top_right[2] > bottom_left[2], "proper rectangle required")
    fe_w = (top_right[1] - bottom_left[1]) / ncols
    fe_h = (top_right[2] - bottom_left[2]) / nrows
    els = nrows * ncols
    nb_vsides = nrows * (ncols-1)
    nb_hsides = (nrows-1) * ncols
    nb_sides = nb_vsides + nb_hsides
    sidenum_first_nb_hside = nb_vsides + 1

    fe_dims = [fe_w, fe_h]
    fe_dims_wo_1 = [fe_h]
    fe_dims_wo_2 = [fe_w]

    new(bottom_left[1], bottom_left[2], top_right[1], top_right[2],
        convert(MeshCoord, nrows), convert(MeshCoord, ncols),
        els,
        nb_sides,
        nb_vsides,
        nb_hsides,
        sidenum_first_nb_hside,
        fe_dims,
        fe_dims_wo_1,
        fe_dims_wo_2,
        Array(R, 2), # integrand args work array
        zeros(R, 2), # ref fe min bounds
        zeros(R, 1), # ref fe min bounds, short
        integration_rel_err,
        integration_abs_err)
  end # RectMesh2 constructor

end # type RectMesh2

RectMesh2(bottom_left::(R,R), top_right::(R,R), nrows::Integer, ncols::Integer) =
  RectMesh2(bottom_left, top_right, nrows, ncols, default_integration_rel_err, default_integration_abs_err)


# ------------------------------------------
# Implement functions required of all meshes.

import Mesh.space_dim
space_dim(m::RectMesh2) = dim(2)

const one_monomial = Monomial(0,0)
import Mesh.one_mon
one_mon(m::RectMesh2) = one_monomial

import Mesh.num_fes
num_fes(m::RectMesh2) = m.num_elements

import Mesh.num_nb_sides
num_nb_sides(m::RectMesh2) = m.num_nb_sides

import Mesh.num_side_faces_per_fe
num_side_faces_per_fe(mesh::RectMesh2) = 4

import Mesh.dependent_dim_for_nb_side
dependent_dim_for_nb_side(i::NBSideNum, mesh::RectMesh2) =
  is_vert_nb_side(i, mesh) ? 1 : 2

import Mesh.dependent_dim_for_ref_side_face
dependent_dim_for_ref_side_face(side_face::FEFace, mesh::RectMesh2) =
  side_face == left_face || side_face == right_face ? 1 : 2


import Mesh.fe_inclusions_of_nb_side!
function fe_inclusions_of_nb_side!(i::NBSideNum, mesh::RectMesh2, fe_incls::NBSideInclusions)
  const mesh_cols = mesh.cols
  if is_vert_nb_side(i, mesh)
    # See "vertical nb sides ordering" diagram above as a reference for the vertical side numbering.
    # Convert to side row and column indexes (only considering the vertical sides, so rows have mesh_cols-1 items).
    const vside_ix = i - 1
    const vside_cols = mesh_cols - 1
    const vside_rix = div(vside_ix, vside_cols)
    const vside_cix = mod(vside_ix, vside_cols)
    const fes_in_complete_rows_prec = vside_rix * mesh_cols # num of fes in complete rows preceding the side's row
    # Find finite element numbers on left and right of the side (remember the fe rows have one more element than the side rows).
    fe_incls.fe1 = fes_in_complete_rows_prec + vside_cix + 1 # the fe on the left of the side
    fe_incls.face_in_fe1 = right_face # on the left fe the side is the right face
    fe_incls.fe2 = fes_in_complete_rows_prec + vside_cix + 2 # the fe on the right of the side
    fe_incls.face_in_fe2 = left_face # on the right fe the side is the left face
  elseif is_horz_nb_side(i, mesh)
    # See the "horizontal nb sides ordering" diagram above as a reference for the non-boundary horizontal side ordering.
    const hside_num = i - mesh.num_nb_vert_sides
    fe_incls.fe1 = hside_num # bottom fe
    fe_incls.face_in_fe1 = top_face # on the bottom fe the side is the top face
    fe_incls.fe2  = hside_num + mesh.cols # top fe
    fe_incls.face_in_fe2 = bottom_face # on the top fe the bel support is the bottom face
  else
    error("side number out of range")
  end
end

import Mesh.is_boundary_side
function is_boundary_side(fe::FENum, face::FEFace, mesh::RectMesh2)
  if face == top_face
    fe_row(fe, mesh) == mesh.rows
  elseif face == bottom_face
    fe_row(fe, mesh) == 1
  elseif face == right_face
    fe_col(fe, mesh) == mesh.cols
  elseif face == left_face
    fe_col(fe, mesh) == 1
  else
    error("invalid face: $face")
  end
end

# TODO: These should use the *side-local* coordinate system in each case below.  In which case it brings up the
#       question of why we care about the side's role ("face") in the reference fe at all?
import Mesh.integral_on_ref_fe_face
function integral_on_ref_fe_face(mon::Monomial, face::FEFace, mesh::RectMesh2)
  if face == Mesh.interior_face
    Poly.integral_on_rect_at_origin(mon, mesh.fe_dims)
  elseif face == top_face
    dim_red_intgd = Poly.reduce_dim_by_fixing(dim(2), mesh.fe_dims[2], mon)
    Poly.integral_on_rect_at_origin(dim_red_intgd, mesh.fe_dims_wo_2) # integral over dim 1 from 0 to fe_width
  elseif face == right_face
    dim_red_intgd = Poly.reduce_dim_by_fixing(dim(1), mesh.fe_dims[1], mon)
    Poly.integral_on_rect_at_origin(dim_red_intgd, mesh.fe_dims_wo_1) # integral over the remaining non-fixed dimension
  elseif face == bottom_face
    if mon.exps[2] != 0
      zeroR
    else
      dim_red_intgd = Poly.reduce_dim_by_fixing(dim(2), zeroR, mon)
      Poly.integral_on_rect_at_origin(dim_red_intgd, mesh.fe_dims_wo_2) # integral over the remaining non-fixed dimension
    end
  elseif face == left_face
    if mon.exps[1] != 0
      zeroR
    else
      dim_red_intgd = Poly.reduce_dim_by_fixing(dim(1), zeroR, mon)
      Poly.integral_on_rect_at_origin(dim_red_intgd, mesh.fe_dims_wo_1) # integral over the remaining non-fixed dimensions
    end
  else
    error("invalid face: $face")
  end
end

import Mesh.integral_on_ref_fe_side_vs_outward_normal
function integral_on_ref_fe_side_vs_outward_normal(vm::VectorMonomial, face::FEFace, mesh::RectMesh2)
  if face == top_face
    integral_on_ref_fe_face(vm[2], face, mesh)
  elseif face == right_face
    integral_on_ref_fe_face(vm[1], face, mesh)
  elseif face == bottom_face
    -integral_on_ref_fe_face(vm[2], face, mesh)
  elseif face == left_face
    -integral_on_ref_fe_face(vm[1], face, mesh)
  else
    error("invalid face: $face")
  end
end

import Mesh.integral_prod_on_fe_face
function integral_prod_on_fe_face(f::Function, mon::Monomial, fe::FENum, face::FEFace, mesh::RectMesh2)
  const d = 2
  const fe_local_origin = fe_coords(fe, mesh)
  const fe_x = mesh.intgd_args_work_array
  if face == Mesh.interior_face
    function ref_intgd(x::Vector{R})
      fe_x[1] = fe_local_origin[1] + x[1]
      fe_x[2] = fe_local_origin[2] + x[2]
      f(fe_x) * Poly.monomial_value(mon, x)
    end
    hcubature(ref_intgd, mesh.ref_fe_min_bounds, mesh.fe_dims, mesh.integration_rel_err, mesh.integration_abs_err)[1]
  else # side face
    # perp axis for the side
    const a = face == left_face || face == right_face ? 1 : 2
    const a_coord_of_ref_side = face == left_face || face == bottom_face ? zeroR : mesh.fe_dims[a]
    const a_coord_of_fe_side = fe_local_origin[a] + a_coord_of_ref_side
    const dim_reduced_poly = Poly.reduce_dim_by_fixing(dim(a), a_coord_of_ref_side, mon)
    function ref_intgd(x::R)
      if a == 1
        fe_x[1] = a_coord_of_fe_side
        fe_x[2] = fe_local_origin[2] + x
      else
        fe_x[1] = fe_local_origin[1] + x
        fe_x[2] = a_coord_of_fe_side
      end
      f(fe_x) * Poly.polynomial_value(dim_reduced_poly, x)
    end
    const non_a_fe_dim = a == 1 ? mesh.fe_dims_wo_1[1] : mesh.fe_dims_wo_2[1]
    hquadrature(ref_intgd, 0., non_a_fe_dim, mesh.integration_rel_err, mesh.integration_abs_err)[1]
  end
end


#
# ------------------------------------------

# functions specific to rectangular meshes

is_vert_nb_side(i::NBSideNum, mesh::RectMesh2) = i < mesh.sidenum_first_nb_horz_side
is_horz_nb_side(i::NBSideNum, mesh::RectMesh2) = mesh.sidenum_first_nb_horz_side <= i <= mesh.num_nb_sides

# face numbers
# Only side faces are defined here, the interior_face is defined in the Mesh module.
const top_face      = fe_face(1)
const right_face    = fe_face(2)
const bottom_face   = fe_face(3)
const left_face     = fe_face(4)
const num_sides = 0x04

fe_row(fe::FENum, mesh::RectMesh2) = fe_rix(fe-1, mesh) + 1
fe_col(fe::FENum, mesh::RectMesh2) = fe_cix(fe-1, mesh) + 1

fe_rix(fe_ix::FENum, mesh::RectMesh2) = div(fe_ix, mesh.cols)
fe_cix(fe_ix::FENum, mesh::RectMesh2) = mod(fe_ix, mesh.cols)

# Fill the passed 2 element array with the coordinates of the corner with range minimum coordinates.
function fe_coords!(fe::FENum, mesh::RectMesh2, coords::Vector{R})
  const fe_ix = fe - 1
  coords[1] = mesh.bottom_left_x + fe_cix(fe_ix, mesh) * mesh.fe_dims[1]
  coords[2] = mesh.bottom_left_y + fe_rix(fe_ix, mesh) * mesh.fe_dims[2]
end

# Functional variant of the above.
fe_coords(fe::FENum, mesh::RectMesh2) =
  let a = zeros(R,2)
    fe_coords!(fe, mesh, a)
    a
  end

end # end of module
