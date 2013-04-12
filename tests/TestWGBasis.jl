using Test
require("../Common")
require("../Poly")
require("../Mesh")
require("../WGrad")
require("../WGBasis")
require("../RMesh")

using WGBasis
using Common
import Mesh, Mesh.fe_num
import RMesh, RMesh.mesh_coord
import Poly, Poly.Monomial
import WGrad, WGrad.WGradSolver

# 2 x 3 mesh, k = 2
#  ----------
#  |  |  |  |
#  ----------
#  |  |  |  |
#  ----------
#
# 6 interior monomials, 2 side monomials
# interior basis els: 3 * 2 * 6 = 36
# vertical side basis els: 2 * 2 * 2 = 8
# horizontal side basis els: 3 * 2 = 6

basis = WeakFunsPolyBasis(deg(2), deg(1), RMesh.RectMesh([0.,0.], [3.,2.], [mesh_coord(3),mesh_coord(2)]))

x = Monomial(1,0)
y = Monomial(0,1)
one = Monomial(0,0)

@test map(m -> m.exps, ref_interior_mons(basis)) == {[0x00, 0x00], [0x00, 0x01], [0x00, 0x02], [0x01, 0x00], [0x01, 0x01], [0x02, 0x00]}

# reference mons for sides for which dimension 1 is dependent on the others
@test basis.ref_side_mons[1] == [one, y]
# reference mons for sides for which dimension 2 is dependent on the others
@test basis.ref_side_mons[2] == [one, x]

@test length(ref_interior_mons(basis)) == 6 == Poly.count_monomials_of_degree_le(deg(2), dim(2))

# sides have reduced domain dimensions, since some dimension is affine-dependent on the others
@test length(basis.ref_side_mons[1]) == 2 == Poly.count_monomials_of_degree_le(deg(1), dim(1))
@test length(basis.ref_side_mons[2]) == 2 == Poly.count_monomials_of_degree_le(deg(1), dim(1))
@test length(basis.ref_side_mons[1]) == basis.mons_per_fe_side

@test length(ref_interior_mons(basis)) == basis.mons_per_fe_interior

@test is_interior_supported(bel_num(1), basis)
@test is_interior_supported(bel_num(36), basis)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(37), basis), basis.mesh) == dim(1)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(44), basis), basis.mesh) == dim(1)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(45), basis), basis.mesh) == dim(2)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(49), basis), basis.mesh) == dim(2)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(50), basis), basis.mesh) == dim(2)
@test_fails RMesh.perp_axis_for_nb_side(support_nb_side_num(bel_num(51), basis), basis.mesh)

function support_fes(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_supported(i, basis)
    [support_interior_num(i, basis)]
  else
    fe_incls = fe_inclusions_of_side_support(i, basis)
    [fe_incls.fe1, fe_incls.fe2]
  end
end

# The first section of basis elements are monomials which are assigned to interiors in blocks of 6.
@test support_fes(bel_num(1), basis) == [1]
@test support_fes(bel_num(6), basis) == [1]
@test support_fes(bel_num(7), basis) == [2]
@test support_fes(bel_num(12), basis) == [2]
@test support_fes(bel_num(13), basis) == [3]
@test support_fes(bel_num(30), basis) == [5]
@test support_fes(bel_num(31), basis) == [6]
@test support_fes(bel_num(36), basis) == [6]

# The next section of basis elements represent monomials on vertical sides, in blocks of 2.
@test support_fes(bel_num(37), basis) == [1,2]
@test support_fes(bel_num(38), basis) == [1,2]
@test support_fes(bel_num(39), basis) == [2,3]
@test support_fes(bel_num(40), basis) == [2,3]
@test support_fes(bel_num(41), basis) == [4,5]
@test support_fes(bel_num(42), basis) == [4,5]
@test support_fes(bel_num(43), basis) == [5,6]
@test support_fes(bel_num(44), basis) == [5,6]

# The final section of basis elements are the horizontal side monomials, assigned to sides in blocks of 3.
@test support_fes(bel_num(45), basis) == [1,4]
@test support_fes(bel_num(46), basis) == [1,4]
@test support_fes(bel_num(47), basis) == [2,5]
@test support_fes(bel_num(48), basis) == [2,5]
@test support_fes(bel_num(49), basis) == [3,6]
@test support_fes(bel_num(50), basis) == [3,6]

left_face = RMesh.lesser_side_face_perp_to_axis(dim(1))
right_face = RMesh.greater_side_face_perp_to_axis(dim(1))
bottom_face = RMesh.lesser_side_face_perp_to_axis(dim(2))
top_face = RMesh.greater_side_face_perp_to_axis(dim(2))

# fe vertical side supported basis element
incls = Mesh.fe_inclusions_of_nb_side(support_nb_side_num(bel_num(37), basis), basis.mesh)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == left_face

incls = fe_inclusions_of_side_support(bel_num(37), basis)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == left_face

# fe horizontal side supported basis element
incls = fe_inclusions_of_side_support(bel_num(45), basis)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == bottom_face


# Test retrieving side monomial number for given monomial and face.
@test mon_num_for_mon_on_side_face(one, right_face, basis) == 1
@test mon_num_for_mon_on_side_face(one, left_face, basis)  == 1
@test mon_num_for_mon_on_side_face(y,   right_face, basis) == 2
@test mon_num_for_mon_on_side_face(y,   left_face, basis)  == 2
@test mon_num_for_mon_on_side_face(one, top_face, basis)    == 1
@test mon_num_for_mon_on_side_face(one, bottom_face, basis) == 1
@test mon_num_for_mon_on_side_face(x,   top_face, basis)    == 2
@test mon_num_for_mon_on_side_face(x,   bottom_face, basis) == 2

# Test support face monomial retrieval

# all interior monomials in finite element 1
@test interior_monomial(bel_num(1), basis) == one
@test interior_monomial(bel_num(2), basis) == y
@test interior_monomial(bel_num(3), basis) == y^2
@test interior_monomial(bel_num(4), basis) == x
@test interior_monomial(bel_num(5), basis) == x*y
@test interior_monomial(bel_num(6), basis) == x^2

@test ref_interior_mons(basis) == [one, y, y^2, x, x*y, x^2]

@test interior_monomial_num(bel_num(3), basis) == mon_num(3)
@test interior_monomial_num(bel_num(6), basis) == mon_num(6)

@test interior_monomial_by_num(mon_num(3), basis) == y^2
@test interior_monomial_by_num(mon_num(6), basis) == x^2


# interior monomials, finite element 2
@test interior_monomial(bel_num(7), basis) == one
@test interior_monomial(bel_num(12), basis) == x^2

@test interior_monomial_num(bel_num(7), basis) == 1
@test interior_monomial_num(bel_num(12), basis) == 6

# interior monomials, finite element 6
@test interior_monomial(bel_num(31), basis) == one
@test interior_monomial(bel_num(36), basis) == x^2

@test interior_monomial_num(bel_num(31), basis) == 1
@test interior_monomial_num(bel_num(36), basis) == 6

# monomials on vertical side between finite elements 1 and 2
@test side_monomial(bel_num(37), basis) == one
@test side_monomial(bel_num(38), basis) == y

@test side_monomial_num(bel_num(37), basis) == 1
@test side_monomial_num(bel_num(38), basis) == 2

@test side_monomial_by_face_and_num(right_face, mon_num(1), basis) == side_monomial_by_face_and_num(left_face, mon_num(1), basis) == one
@test side_monomial_by_face_and_num(right_face, mon_num(2), basis) == side_monomial_by_face_and_num(right_face, mon_num(2), basis) == y
@test ref_side_mons(right_face, basis) == ref_side_mons(left_face, basis) == [one, y]
@test_fails side_monomial_by_face_and_num(right_face, mon_num(3), basis)

# monomials on vertical side between finite elements 2 and 3
@test side_monomial(bel_num(39), basis) == one
@test side_monomial(bel_num(40), basis) == y

# monomials on vertical side between finite elements 5 and 6
@test side_monomial(bel_num(43), basis) == one
@test side_monomial(bel_num(44), basis) == y

# monomials on horizontal side between finite elements 1 and 4
@test side_monomial(bel_num(45), basis) == one
@test side_monomial(bel_num(46), basis) == x

@test side_monomial_by_face_and_num(top_face, mon_num(1), basis) == side_monomial_by_face_and_num(bottom_face, mon_num(1), basis) == one
@test side_monomial_by_face_and_num(top_face, mon_num(2), basis) == side_monomial_by_face_and_num(bottom_face, mon_num(2), basis) == x
@test ref_side_mons(top_face, basis) == ref_side_mons(bottom_face, basis) == [one, x]
@test_fails side_monomial_by_face_and_num(bottom_face, mon_num(3), basis)

# monomials on horizontal side between finite elements 2 and 5
@test side_monomial(bel_num(47), basis) == one
@test side_monomial(bel_num(48), basis) == x

# monomials on horizontal side between finite elements 3 and 6
@test side_monomial(bel_num(49), basis) == one
@test side_monomial(bel_num(50), basis) == x


# test weak gradients of basis elements

wgrad_solver = WGradSolver(deg(1), basis.mesh)

# weak gradients of interior supported elements
@test wgrad_for_mon_on_interior(mon_num(1), basis) == WGrad.wgrad(one, Mesh.interior_face, wgrad_solver)
@test wgrad_for_mon_on_interior(mon_num(2), basis) == WGrad.wgrad(y, Mesh.interior_face, wgrad_solver)
@test wgrad_for_mon_on_interior(mon_num(3), basis) == WGrad.wgrad(y^2, Mesh.interior_face, wgrad_solver)
@test wgrad_for_mon_on_interior(mon_num(4), basis) == WGrad.wgrad(x, Mesh.interior_face, wgrad_solver)
@test wgrad_for_mon_on_interior(mon_num(5), basis) == WGrad.wgrad(x*y, Mesh.interior_face, wgrad_solver)
@test wgrad_for_mon_on_interior(mon_num(6), basis) == WGrad.wgrad(x^2, Mesh.interior_face, wgrad_solver)

@test wgrad_for_mon_on_side_face(mon_num(1), right_face, basis) == WGrad.wgrad(one, right_face, wgrad_solver)
@test wgrad_for_mon_on_side_face(mon_num(2), right_face, basis) == WGrad.wgrad(y, right_face, wgrad_solver)
@test_fails wgrad_for_mon_on_side_face(mon_num(3), right_face, basis)

# weak gradients for monomials on vertical side between finite elements 1 and 2
@test wgrads_for_side_bel(bel_num(37), basis).wgrad_on_fe1 == WGrad.wgrad(one, right_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(37), basis).wgrad_on_fe2 == WGrad.wgrad(one, left_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(38), basis).wgrad_on_fe1 == WGrad.wgrad(y, right_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(38), basis).wgrad_on_fe2 == WGrad.wgrad(y, left_face, wgrad_solver)


# weak gradients for monomials on horizontal side between finite elements 1 and 4
@test wgrads_for_side_bel(bel_num(45), basis).wgrad_on_fe1 == WGrad.wgrad(one, top_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(45), basis).wgrad_on_fe2 == WGrad.wgrad(one, bottom_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(46), basis).wgrad_on_fe1 == WGrad.wgrad(x, top_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(46), basis).wgrad_on_fe2 == WGrad.wgrad(x, bottom_face, wgrad_solver)


# test L2 inner products of basis elements

# interior supported element inner products
int_ips = ips_ref_interior_mons(basis)
@test int_ips[1,1] == Mesh.integral_face_rel_x_face_rel_on_face(one, one, Mesh.interior_face, basis.mesh) == 1
@test int_ips[1,2] == Mesh.integral_face_rel_x_face_rel_on_face(one, y, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,3] == Mesh.integral_face_rel_x_face_rel_on_face(one, y^2, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[1,4] == Mesh.integral_face_rel_x_face_rel_on_face(one, x, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,5] == Mesh.integral_face_rel_x_face_rel_on_face(one, x*y, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[1,6] == Mesh.integral_face_rel_x_face_rel_on_face(one, x^2, Mesh.interior_face, basis.mesh) == 1/3

@test int_ips[2,1] == Mesh.integral_face_rel_x_face_rel_on_face(y, one, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[2,2] == Mesh.integral_face_rel_x_face_rel_on_face(y, y, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[2,3] == Mesh.integral_face_rel_x_face_rel_on_face(y, y^2, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,4] == Mesh.integral_face_rel_x_face_rel_on_face(y, x, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,5] == Mesh.integral_face_rel_x_face_rel_on_face(y, x*y, Mesh.interior_face, basis.mesh) == 1/6
@test int_ips[2,6] == Mesh.integral_face_rel_x_face_rel_on_face(y, x^2, Mesh.interior_face, basis.mesh) == 1/6

# side supported element inner products
top_face_ips = ips_ref_side_mons(top_face, basis)

@test top_face_ips[1,1] == Mesh.integral_face_rel_x_face_rel_on_face(one, one, top_face, basis.mesh) == 1
@test top_face_ips[1,2] == Mesh.integral_face_rel_x_face_rel_on_face(one, x, top_face, basis.mesh) == 1/2
@test top_face_ips[2,1] == top_face_ips[1,2]
@test top_face_ips[2,2] == Mesh.integral_face_rel_x_face_rel_on_face(x, x, top_face, basis.mesh) == 1/3
