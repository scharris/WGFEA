using Test
require("../Common")
require("../Poly")
require("../Mesh")
require("../RMesh2")

import RMesh2, RMesh2.RectMesh2, RMesh2.fe_col, RMesh2.fe_row, RMesh2.fe_coords
using Common
using Mesh
import Poly.Monomial, Poly.VectorMonomial

rmesh200x400 = RectMesh2((0.,0.), (200.,100.), 200, 400)

# Returns the number of the finite element to the left, or 0 if there is none.
fe_left_of(fe::FENum, mesh::RectMesh2) =
  let col = fe_col(fe, mesh)
    if col == 1 0 else fe - 1 end
  end

# Returns the number of the finite element to the right, or 0 if there is none.
fe_right_of(fe::FENum, mesh::RectMesh2) =
  let col = fe_col(fe, mesh)
    if col == mesh.cols 0 else fe + 1 end
  end

# Returns the number of the finite element above the passed finite element, or 0 if there is none.
fe_above(fe::FENum, mesh::RectMesh2) =
  let row = fe_row(fe, mesh)
    if row == mesh.rows 0 else fe + mesh.cols end
  end

# Returns the number of the finite element below the passed finite element, or 0 if there is none.
fe_below(fe::FENum, mesh::RectMesh2) =
  let row = fe_row(fe, mesh)
    if row == 1 0 else fe - mesh.cols end
  end

@test fe_left_of(fe_num(1), rmesh200x400) == 0
@test fe_left_of(fe_num(2), rmesh200x400) == 1
@test fe_left_of(fe_num(400), rmesh200x400) == 399
@test fe_left_of(fe_num(401), rmesh200x400) == 0
@test fe_left_of(fe_num(199*400+1), rmesh200x400) == 0
@test fe_left_of(fe_num(200*400), rmesh200x400) == 200*400-1

@test fe_right_of(fe_num(1), rmesh200x400) == 2
@test fe_right_of(fe_num(400), rmesh200x400) == 0
@test fe_right_of(fe_num(401), rmesh200x400) == 402
@test fe_right_of(fe_num(199*400+1), rmesh200x400) == 199*400+2
@test fe_right_of(fe_num(200*400), rmesh200x400) == 0

@test fe_above(fe_num(1), rmesh200x400) == 401
@test fe_above(fe_num(401), rmesh200x400) == 801
@test fe_above(fe_num(198*400+1), rmesh200x400) == 199*400+1
@test fe_above(fe_num(200*400), rmesh200x400) == 0

@test fe_below(fe_num(401), rmesh200x400) == 1
@test fe_below(fe_num(400), rmesh200x400) == 0
@test fe_below(fe_num(801), rmesh200x400) == 401
@test fe_below(fe_num(199*400+1), rmesh200x400) == 198*400+1
@test fe_below(fe_num(1), rmesh200x400) == 0


# Test finite element plane coordinates

# 10 rows x 20 cols mesh, with vertexes at each integer pair
rmesh10x20 = RectMesh2((0.,0.), (20.,10.), 10, 20)

@test fe_coords(fe_num(1), rmesh10x20) == [0., 0.]
@test fe_coords(fe_num(2), rmesh10x20) == [1., 0.]
@test fe_coords(fe_num(20), rmesh10x20) == [19., 0.]
@test fe_coords(fe_num(21), rmesh10x20) == [0., 1.]
@test fe_coords(fe_num(200), rmesh10x20) == [19., 9.]

@test rmesh10x20.bottom_left_x == 0.
@test rmesh10x20.bottom_left_y == 0.
@test rmesh10x20.top_right_x == 20.
@test rmesh10x20.top_right_y == 10.
@test rmesh10x20.rows == 10
@test rmesh10x20.cols == 20
@test rmesh10x20.fe_dims[1] == 1.
@test rmesh10x20.fe_dims[2] == 1.
@test rmesh10x20.num_elements == 200
@test rmesh10x20.num_nb_vert_sides == 190
@test rmesh10x20.num_nb_horz_sides == 180
@test rmesh10x20.num_nb_sides == 370
@test rmesh10x20.sidenum_first_nb_horz_side == 191

@test num_fes(rmesh10x20) == 200
@test num_nb_sides(rmesh10x20) == 370

# test boundary sides

@test !is_boundary_side(fe_num(1), RMesh2.top_face, rmesh10x20)
@test !is_boundary_side(fe_num(1), RMesh2.right_face, rmesh10x20)
@test is_boundary_side(fe_num(1), RMesh2.bottom_face, rmesh10x20)
@test is_boundary_side(fe_num(1), RMesh2.left_face, rmesh10x20)

@test !is_boundary_side(fe_num(8), RMesh2.top_face, rmesh10x20)
@test !is_boundary_side(fe_num(8), RMesh2.right_face, rmesh10x20)
@test is_boundary_side(fe_num(8), RMesh2.bottom_face, rmesh10x20)
@test !is_boundary_side(fe_num(8), RMesh2.left_face, rmesh10x20)

@test !is_boundary_side(fe_num(20), RMesh2.top_face, rmesh10x20)
@test is_boundary_side(fe_num(20), RMesh2.right_face, rmesh10x20)
@test is_boundary_side(fe_num(20), RMesh2.bottom_face, rmesh10x20)
@test !is_boundary_side(fe_num(20), RMesh2.left_face, rmesh10x20)

@test !is_boundary_side(fe_num(21), RMesh2.top_face, rmesh10x20)
@test !is_boundary_side(fe_num(21), RMesh2.right_face, rmesh10x20)
@test !is_boundary_side(fe_num(21), RMesh2.bottom_face, rmesh10x20)
@test is_boundary_side(fe_num(21), RMesh2.left_face, rmesh10x20)

@test !is_boundary_side(fe_num(22), RMesh2.top_face, rmesh10x20)
@test !is_boundary_side(fe_num(22), RMesh2.right_face, rmesh10x20)
@test !is_boundary_side(fe_num(22), RMesh2.bottom_face, rmesh10x20)
@test !is_boundary_side(fe_num(22), RMesh2.left_face, rmesh10x20)

@test is_boundary_side(fe_num(200), RMesh2.top_face, rmesh10x20)
@test is_boundary_side(fe_num(200), RMesh2.right_face, rmesh10x20)
@test !is_boundary_side(fe_num(200), RMesh2.bottom_face, rmesh10x20)
@test !is_boundary_side(fe_num(200), RMesh2.left_face, rmesh10x20)


rmesh2x3 = RMesh2.RectMesh2((0.,0.), (3.,2.), 2, 3)
@test RMesh2.is_vert_nb_side(nb_side_num(1), rmesh2x3)
@test RMesh2.is_vert_nb_side(nb_side_num(4), rmesh2x3)
@test RMesh2.is_horz_nb_side(nb_side_num(5), rmesh2x3)
@test RMesh2.is_horz_nb_side(nb_side_num(7), rmesh2x3)
@test RMesh2.is_horz_nb_side(nb_side_num(7), rmesh2x3)
@test !RMesh2.is_horz_nb_side(nb_side_num(8), rmesh2x3)

# Test side inclusions
incls = Mesh.NBSideInclusions()

# fe vertical sides
fe_inclusions_of_nb_side!(nb_side_num(1), rmesh2x3, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(nb_side_num(2), rmesh2x3, incls)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(3)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(nb_side_num(3), rmesh2x3, incls)
@test incls.fe1 == fe_num(4)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(nb_side_num(4), rmesh2x3, incls)
@test incls.fe1 == fe_num(5)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == RMesh2.left_face

# fe horizontal sides
fe_inclusions_of_nb_side!(nb_side_num(5), rmesh2x3, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == RMesh2.bottom_face

fe_inclusions_of_nb_side!(nb_side_num(6), rmesh2x3, incls)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == RMesh2.bottom_face

fe_inclusions_of_nb_side!(nb_side_num(7), rmesh2x3, incls)
@test incls.fe1 == fe_num(3)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == RMesh2.bottom_face



# Test integrals on finite element faces.

x = Monomial(1,0)
y = Monomial(0,1)
one_mon = RMesh2.one_mon(rmesh2x3)

# Integrate constants on faces of reference element.
@test Mesh.integral_face_rel_on_face(one_mon, Mesh.interior_face, rmesh2x3) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, RMesh2.left_face, rmesh2x3) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, RMesh2.right_face, rmesh2x3) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, RMesh2.bottom_face, rmesh2x3) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, RMesh2.top_face, rmesh2x3) == 1.
@test Mesh.integral_face_rel_on_face(3., Mesh.interior_face, rmesh2x3) == 3.
@test Mesh.integral_face_rel_on_face(3., RMesh2.top_face, rmesh2x3) == 3.

# Integrate monomials on faces of reference element.

@test Mesh.integral_face_rel_on_face(x*y^2, Mesh.interior_face, rmesh2x3) == 1/6
@test Mesh.integral_face_rel_on_face(x*y^2, RMesh2.left_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(x*y^2, RMesh2.right_face, rmesh2x3) == 1/3
@test Mesh.integral_face_rel_on_face(x*y^2, RMesh2.bottom_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(x*y^2, RMesh2.top_face, rmesh2x3) == 1/2

@test Mesh.integral_face_rel_on_face(x^3*y^4, Mesh.interior_face, rmesh2x3) == 1/20
@test Mesh.integral_face_rel_on_face(x^3*y^4, RMesh2.left_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(x^3*y^4, RMesh2.right_face, rmesh2x3) == 1/5
@test Mesh.integral_face_rel_on_face(x^3*y^4, RMesh2.bottom_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(x^3*y^4, RMesh2.top_face, rmesh2x3) == 1/4

## Integrate polynomial on faces of reference element.
@test Mesh.integral_face_rel_on_face(2*x*y^2 + 3*x^3*y^4, Mesh.interior_face, rmesh2x3) == 2*(1/6) + 3*(1/20)
@test Mesh.integral_face_rel_on_face(1.2*x*y^2 + 3.4*x^3*y^4, RMesh2.left_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(1.5*x*y^2 + 0.123*x^3*y^4, RMesh2.right_face, rmesh2x3) == 1.5*(1/3) + 0.123*(1/5)
@test Mesh.integral_face_rel_on_face(x*y^2 + x^3*y^4, RMesh2.bottom_face, rmesh2x3) == 0
@test Mesh.integral_face_rel_on_face(4.5*x*y^2 + x^3*y^4, RMesh2.top_face, rmesh2x3) == 4.5*(1/2) + 1/4

vm1 = VectorMonomial(x^3*y^4, dim(1))
vm2 = VectorMonomial(x^3*y^4, dim(2))

# Integrate vector monomials vs outward normals on the side faces of the reference finite element.

@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh2.left_face, rmesh2x3) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh2.right_face, rmesh2x3) == 1/5
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh2.bottom_face, rmesh2x3) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh2.top_face, rmesh2x3) == 0

@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh2.left_face, rmesh2x3) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh2.right_face, rmesh2x3) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh2.bottom_face, rmesh2x3) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh2.top_face, rmesh2x3) == 1/4

# Integrate on lesser faces using monomials which are constant in the corresponding dimension (so the integrals aren't 0).
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(VectorMonomial(y^4,dim(1)), RMesh2.left_face, rmesh2x3) == -1/5
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(VectorMonomial(x^3,dim(2)), RMesh2.bottom_face, rmesh2x3) == -1/4

nearly_eq(a::R, b::R) = abs(b - a) < 10e-8

# Test integration of a product of an arbitrary function and an element-local monomial on finite element faces.

# Integrating the product below on fe 1 should be equivalent to integrating the monomial x^3 y^4 on the reference element interior.
f1(x::Vector{R}) = x[1]^2 * x[2]^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), Mesh.interior_face, rmesh2x3), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), RMesh2.left_face, rmesh2x3), 0.)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), RMesh2.right_face, rmesh2x3), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), RMesh2.bottom_face, rmesh2x3), 0.)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), RMesh2.top_face, rmesh2x3), 1/4)

# Integrating the product below on fe 5 should be equivalent to integrating the monomial x^3 y^4 z on the reference element interior.
fe5_coords = RMesh2.fe_coords(fe_num(5), rmesh2x3)
f2(x::Vector{R}) = (x[1] - fe5_coords[1])^2 * (x[2] - fe5_coords[2])^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(5), Mesh.interior_face, rmesh2x3), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(5), RMesh2.left_face, rmesh2x3), 0.)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(5), RMesh2.right_face, rmesh2x3), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(5), RMesh2.bottom_face, rmesh2x3), 0.)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(5), RMesh2.top_face, rmesh2x3), 1/4)

