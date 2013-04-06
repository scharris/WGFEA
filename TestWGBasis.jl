using Test
using WGBasis

using Common
import Mesh, Mesh.fe_num
import RMesh2
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

basis = WeakFunsPolyBasis(deg(2), deg(1), RMesh2.RectMesh2((0.,0.), (3.,2.), 2, 3))

x = Monomial(1,0)
y = Monomial(0,1)
one = Monomial(0,0)

@test map(m -> m.exps, basis.ref_interior_mons) == {[0x00, 0x00], [0x00, 0x01], [0x00, 0x02], [0x01, 0x00], [0x01, 0x01], [0x02, 0x00]}

# reference mons for sides for which dimension 1 is dependent on the others
@test basis.ref_side_mons[1] == [one, y]
# reference mons for sides for which dimension 2 is dependent on the others
@test basis.ref_side_mons[2] == [one, x]

@test length(basis.ref_interior_mons) == 6 == Poly.count_monomials_of_degree_le(deg(2), dim(2))

# sides have reduced domain dimensions, since some dimension is affine-dependent on the others
@test length(basis.ref_side_mons[1]) == 2 == Poly.count_monomials_of_degree_le(deg(1), dim(1))
@test length(basis.ref_side_mons[2]) == 2 == Poly.count_monomials_of_degree_le(deg(1), dim(1))
@test length(basis.ref_side_mons[1]) == basis.mons_per_fe_side

@test length(basis.ref_interior_mons) == basis.mons_per_fe_interior

@test is_interior_supported(bel_num(1), basis)
@test is_interior_supported(bel_num(36), basis)
@test RMesh2.is_vert_nb_side(support_nb_side_num(bel_num(37), basis), basis.mesh)
@test RMesh2.is_vert_nb_side(support_nb_side_num(bel_num(44), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(45), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(49), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(50), basis), basis.mesh)
@test_fails !RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(51), basis), basis.mesh)


function support_fes(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_supported(i, basis)
    [support_interior_num(i, basis)]
  else
    side_num = support_nb_side_num(i, basis)
    fe_incls = Mesh.fe_inclusions_of_nb_side(side_num, basis.mesh)
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

# Test side inclusions
incls = Mesh.NBSideInclusions()

# fe vertical side supported basis element
Mesh.fe_inclusions_of_nb_side!(support_nb_side_num(bel_num(37), basis), basis.mesh, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == RMesh2.left_face

# fe horizontal side supported basis element
Mesh.fe_inclusions_of_nb_side!(support_nb_side_num(bel_num(45), basis), basis.mesh, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == RMesh2.bottom_face

# Test support face monomial retrieval

# all interior monomials in finite element 1
@test interior_monomial(bel_num(1), basis) == one
@test interior_monomial(bel_num(2), basis) == y
@test interior_monomial(bel_num(3), basis) == y^2
@test interior_monomial(bel_num(4), basis) == x
@test interior_monomial(bel_num(5), basis) == x*y
@test interior_monomial(bel_num(6), basis) == x^2

# interior monomials, finite element 2
@test interior_monomial(bel_num(7), basis) == one
@test interior_monomial(bel_num(12), basis) == x^2

# interior monomials, finite element 6
@test interior_monomial(bel_num(31), basis) == one
@test interior_monomial(bel_num(36), basis) == x^2

# monomials on vertical side between finite elements 1 and 2
@test side_monomial(bel_num(37), basis) == one
@test side_monomial(bel_num(38), basis) == y

# monomials on vertical side between finite elements 2 and 3
@test side_monomial(bel_num(39), basis) == one
@test side_monomial(bel_num(40), basis) == y

# monomials on vertical side between finite elements 5 and 6
@test side_monomial(bel_num(43), basis) == one
@test side_monomial(bel_num(44), basis) == y

# monomials on horizontal side between finite elements 1 and 4
@test side_monomial(bel_num(45), basis) == one
@test side_monomial(bel_num(46), basis) == x

# monomials on horizontal side between finite elements 2 and 5
@test side_monomial(bel_num(47), basis) == one
@test side_monomial(bel_num(48), basis) == x

# monomials on horizontal side between finite elements 3 and 6
@test side_monomial(bel_num(49), basis) == one
@test side_monomial(bel_num(50), basis) == x


# test weak gradients of basis elements

wgrad_solver = WGradSolver(deg(1), basis.mesh)

# weak gradients of interior supported elements
@test wgrad_for_interior_mon_num(mon_num(1), basis) == WGrad.wgrad(one, Mesh.interior_face, wgrad_solver)
@test wgrad_for_interior_mon_num(mon_num(2), basis) == WGrad.wgrad(y, Mesh.interior_face, wgrad_solver)
@test wgrad_for_interior_mon_num(mon_num(3), basis) == WGrad.wgrad(y^2, Mesh.interior_face, wgrad_solver)
@test wgrad_for_interior_mon_num(mon_num(4), basis) == WGrad.wgrad(x, Mesh.interior_face, wgrad_solver)
@test wgrad_for_interior_mon_num(mon_num(5), basis) == WGrad.wgrad(x*y, Mesh.interior_face, wgrad_solver)
@test wgrad_for_interior_mon_num(mon_num(6), basis) == WGrad.wgrad(x^2, Mesh.interior_face, wgrad_solver)

# weak gradients for monomials on vertical side between finite elements 1 and 2
@test wgrads_for_side_bel(bel_num(37), basis).wgrad_on_fe1 == WGrad.wgrad(one, RMesh2.right_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(37), basis).wgrad_on_fe2 == WGrad.wgrad(one, RMesh2.left_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(38), basis).wgrad_on_fe1 == WGrad.wgrad(y, RMesh2.right_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(38), basis).wgrad_on_fe2 == WGrad.wgrad(y, RMesh2.left_face, wgrad_solver)

# weak gradients for monomials on horizontal side between finite elements 1 and 4
@test wgrads_for_side_bel(bel_num(45), basis).wgrad_on_fe1 == WGrad.wgrad(one, RMesh2.top_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(45), basis).wgrad_on_fe2 == WGrad.wgrad(one, RMesh2.bottom_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(46), basis).wgrad_on_fe1 == WGrad.wgrad(x, RMesh2.top_face, wgrad_solver)
@test wgrads_for_side_bel(bel_num(46), basis).wgrad_on_fe2 == WGrad.wgrad(x, RMesh2.bottom_face, wgrad_solver)


# test L2 inner products of basis elements

# interior supported element inner products
int_ips = ips_ref_interior_mons(basis)
@test int_ips[1,1] == Mesh.integral_prod_on_ref_fe_face(one, one, Mesh.interior_face, basis.mesh) == 1
@test int_ips[1,2] == Mesh.integral_prod_on_ref_fe_face(one, y, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,3] == Mesh.integral_prod_on_ref_fe_face(one, y^2, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[1,4] == Mesh.integral_prod_on_ref_fe_face(one, x, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,5] == Mesh.integral_prod_on_ref_fe_face(one, x*y, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[1,6] == Mesh.integral_prod_on_ref_fe_face(one, x^2, Mesh.interior_face, basis.mesh) == 1/3

@test int_ips[2,1] == Mesh.integral_prod_on_ref_fe_face(y, one, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[2,2] == Mesh.integral_prod_on_ref_fe_face(y, y, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[2,3] == Mesh.integral_prod_on_ref_fe_face(y, y^2, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,4] == Mesh.integral_prod_on_ref_fe_face(y, x, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,5] == Mesh.integral_prod_on_ref_fe_face(y, x*y, Mesh.interior_face, basis.mesh) == 1/6
@test int_ips[2,6] == Mesh.integral_prod_on_ref_fe_face(y, x^2, Mesh.interior_face, basis.mesh) == 1/6

# side supported element inner products
top_face_ips = ips_ref_side_mons(RMesh2.top_face, basis)

@test top_face_ips[1,1] == Mesh.integral_prod_on_ref_fe_face(one, one, RMesh2.top_face, basis.mesh) == 1
@test top_face_ips[1,2] == Mesh.integral_prod_on_ref_fe_face(one, x, RMesh2.top_face, basis.mesh) == 1/2
@test top_face_ips[2,1] == top_face_ips[1,2]
@test top_face_ips[2,2] == Mesh.integral_prod_on_ref_fe_face(x, x, RMesh2.top_face, basis.mesh) == 1/3

Q: Are side basis elements well defined?  We are using local monomials defined by the fe's coordinate system,
   but the side has two such fe's, leading to contradictory definitions.
