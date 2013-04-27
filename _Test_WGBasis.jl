using Test
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

# test equality and hashing
@test basis == basis
@test hash(basis) == hash(basis)
basis2 = WeakFunsPolyBasis(deg(2), deg(1), RMesh.RectMesh([0.,0.], [3.,2.], [mesh_coord(3),mesh_coord(2)]))
@test basis == basis2
@test hash(basis) == hash(basis2)
basis2 = WeakFunsPolyBasis(deg(3), deg(1), RMesh.RectMesh([0.,0.], [3.,2.], [mesh_coord(3),mesh_coord(2)]))
@test basis != basis2
basis2 = WeakFunsPolyBasis(deg(2), deg(2), RMesh.RectMesh([0.,0.], [3.,2.], [mesh_coord(3),mesh_coord(2)]))
@test basis != basis2
basis2 = WeakFunsPolyBasis(deg(2), deg(1), RMesh.RectMesh([1.,0.], [3.,2.], [mesh_coord(3),mesh_coord(2)]))
@test basis != basis2
basis2 = WeakFunsPolyBasis(deg(2), deg(1), RMesh.RectMesh([0.,0.], [3.,2.], [mesh_coord(4),mesh_coord(2)]))
@test basis != basis2


x = Monomial(1,0)
y = Monomial(0,1)
one = Monomial(0,0)

@test map(m -> m.exps, interior_mons(basis)) == {[0x00, 0x00], [0x00, 0x01], [0x00, 0x02], [0x01, 0x00], [0x01, 0x01], [0x02, 0x00]}

# reference mons for sides for which dimension 1 is dependent on the others
@test basis.side_mons_by_dep_dim[1] == [one, y]
# reference mons for sides for which dimension 2 is dependent on the others
@test basis.side_mons_by_dep_dim[2] == [one, x]

@test length(interior_mons(basis)) == 6 == Poly.count_mons_of_deg_le(deg(2), dim(2))

# sides have reduced domain dimensions, since some dimension is affine-dependent on the others
@test length(basis.side_mons_by_dep_dim[1]) == 2 == Poly.count_mons_of_deg_le(deg(1), dim(1))
@test length(basis.side_mons_by_dep_dim[2]) == 2 == Poly.count_mons_of_deg_le(deg(1), dim(1))
@test length(basis.side_mons_by_dep_dim[1]) == basis.mons_per_fe_side

@test length(interior_mons(basis)) == basis.mons_per_fe_interior

@test is_interior_supported(beln(1), basis)
@test is_interior_supported(beln(36), basis)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(37), basis), basis.mesh) == dim(1)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(44), basis), basis.mesh) == dim(1)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(45), basis), basis.mesh) == dim(2)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(49), basis), basis.mesh) == dim(2)
@test RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(50), basis), basis.mesh) == dim(2)
@test_fails RMesh.perp_axis_for_nb_side(support_nb_side_num(beln(51), basis), basis.mesh)

function support_fes(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_supported(i, basis)
    [support_interior_num(i, basis)]
  else
    fe_incls = fe_inclusions_of_side_support(i, basis)
    [fe_incls.fe1, fe_incls.fe2]
  end
end

# The first section of basis elements are monomials which are assigned to interiors in blocks of 6.
@test support_fes(beln(1), basis) == [1]
@test support_fes(beln(6), basis) == [1]
@test support_fes(beln(7), basis) == [2]
@test support_fes(beln(12), basis) == [2]
@test support_fes(beln(13), basis) == [3]
@test support_fes(beln(30), basis) == [5]
@test support_fes(beln(31), basis) == [6]
@test support_fes(beln(36), basis) == [6]

# The next section of basis elements represent monomials on vertical sides, in blocks of 2.
@test support_fes(beln(37), basis) == [1,2]
@test support_fes(beln(38), basis) == [1,2]
@test support_fes(beln(39), basis) == [2,3]
@test support_fes(beln(40), basis) == [2,3]
@test support_fes(beln(41), basis) == [4,5]
@test support_fes(beln(42), basis) == [4,5]
@test support_fes(beln(43), basis) == [5,6]
@test support_fes(beln(44), basis) == [5,6]

# The final section of basis elements are the horizontal side monomials, assigned to sides in blocks of 3.
@test support_fes(beln(45), basis) == [1,4]
@test support_fes(beln(46), basis) == [1,4]
@test support_fes(beln(47), basis) == [2,5]
@test support_fes(beln(48), basis) == [2,5]
@test support_fes(beln(49), basis) == [3,6]
@test support_fes(beln(50), basis) == [3,6]

left_face = RMesh.lesser_side_face_perp_to_axis(dim(1))
right_face = RMesh.greater_side_face_perp_to_axis(dim(1))
bottom_face = RMesh.lesser_side_face_perp_to_axis(dim(2))
top_face = RMesh.greater_side_face_perp_to_axis(dim(2))

# fe vertical side supported basis element
incls = Mesh.fe_inclusions_of_nb_side(support_nb_side_num(beln(37), basis), basis.mesh)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == left_face

incls = fe_inclusions_of_side_support(beln(37), basis)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == left_face

# fe horizontal side supported basis element
incls = fe_inclusions_of_side_support(beln(45), basis)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == bottom_face

rshape = Mesh.oshape(1)


# Test support face monomial retrieval

# all interior monomials in finite element 1
@test interior_mon(beln(1), basis) == one
@test interior_mon(beln(2), basis) == y
@test interior_mon(beln(3), basis) == y^2
@test interior_mon(beln(4), basis) == x
@test interior_mon(beln(5), basis) == x*y
@test interior_mon(beln(6), basis) == x^2

@test interior_mons(basis) == [one, y, y^2, x, x*y, x^2]

@test WGBasis.first_bel_supported_on_fe_interior(fe_num(1), basis) == 1

@test interior_mon_num(beln(3), basis) == mon_num(3)
@test interior_mon_num(beln(6), basis) == mon_num(6)

@test interior_mons(basis)[3] == y^2
@test interior_mons(basis)[6] == x^2

# interior monomials, finite element 2
@test WGBasis.first_bel_supported_on_fe_interior(fe_num(2), basis) == 7

@test interior_mon(beln(7), basis) == one
@test interior_mon(beln(12), basis) == x^2

@test interior_mon_num(beln(7), basis) == 1
@test interior_mon_num(beln(12), basis) == 6


# interior monomials, finite element 6
@test WGBasis.first_bel_supported_on_fe_interior(fe_num(6), basis) == 31

@test interior_mon(beln(31), basis) == one
@test interior_mon(beln(36), basis) == x^2

@test interior_mon_num(beln(31), basis) == 1
@test interior_mon_num(beln(36), basis) == 6

# monomials on vertical side between finite elements 1 and 2
@test WGBasis.first_bel_supported_on_fe_side(fe_num(1), right_face, basis) == 37
@test WGBasis.first_bel_supported_on_fe_side(fe_num(2), left_face, basis) == 37

@test side_mon(beln(37), basis) == one
@test side_mon(beln(38), basis) == y

@test side_mon_num(beln(37), basis) == 1
@test side_mon_num(beln(38), basis) == 2

@test side_mons_for_oshape_side(rshape, right_face, basis)[1] == side_mons_for_oshape_side(rshape, left_face, basis)[1] == one
@test side_mons_for_oshape_side(rshape, right_face, basis)[2]  == side_mons_for_oshape_side(rshape, left_face, basis)[2] == y
@test side_mons_for_fe_side(fe_num(1), right_face, basis) == side_mons_for_fe_side(fe_num(1), left_face, basis) == [one, y]
@test_fails side_mons_for_oshape_side(rshape, right_face, basis)[3]

# monomials on vertical side between finite elements 2 and 3
@test WGBasis.first_bel_supported_on_fe_side(fe_num(2), right_face, basis) == 39
@test WGBasis.first_bel_supported_on_fe_side(fe_num(3), left_face, basis) == 39

@test side_mon(beln(39), basis) == one
@test side_mon(beln(40), basis) == y

# monomials on vertical side between finite elements 5 and 6
@test WGBasis.first_bel_supported_on_fe_side(fe_num(5), right_face, basis) == 43
@test WGBasis.first_bel_supported_on_fe_side(fe_num(6), left_face, basis) == 43
@test side_mon(beln(43), basis) == one
@test side_mon(beln(44), basis) == y

# monomials on horizontal side between finite elements 1 and 4
@test WGBasis.first_bel_supported_on_fe_side(fe_num(1), top_face, basis) == 45
@test WGBasis.first_bel_supported_on_fe_side(fe_num(4), bottom_face, basis) == 45
@test side_mon(beln(45), basis) == one
@test side_mon(beln(46), basis) == x

@test side_mons_for_oshape_side(rshape, top_face, basis)[1] == side_mons_for_oshape_side(rshape, bottom_face, basis)[1] == one
@test side_mons_for_oshape_side(rshape, top_face, basis)[2] == side_mons_for_oshape_side(rshape, bottom_face, basis)[2] == x
@test side_mons_for_fe_side(fe_num(1), top_face, basis) == side_mons_for_fe_side(fe_num(1), bottom_face, basis) == [one, x]
@test_fails side_mons_for_oshape_side(rshape, bottom_face, basis)[3]

# monomials on horizontal side between finite elements 2 and 5
@test side_mon(beln(47), basis) == one
@test side_mon(beln(48), basis) == x

# monomials on horizontal side between finite elements 3 and 6
@test side_mon(beln(49), basis) == one
@test side_mon(beln(50), basis) == x


# test weak gradients of basis elements

wgrad_solver = WGradSolver(deg(1), basis.mesh)

# weak gradients of interior supported elements
@test wgrad_interior_mon(mon_num(1), rshape, basis) == WGrad.wgrad(one, rshape, Mesh.interior_face, wgrad_solver)
@test wgrad_interior_mon(mon_num(2), rshape, basis) == WGrad.wgrad(y, rshape, Mesh.interior_face, wgrad_solver)
@test wgrad_interior_mon(mon_num(3), rshape, basis) == WGrad.wgrad(y^2, rshape, Mesh.interior_face, wgrad_solver)
@test wgrad_interior_mon(mon_num(4), rshape, basis) == WGrad.wgrad(x, rshape, Mesh.interior_face, wgrad_solver)
@test wgrad_interior_mon(mon_num(5), rshape, basis) == WGrad.wgrad(x*y, rshape, Mesh.interior_face, wgrad_solver)
@test wgrad_interior_mon(mon_num(6), rshape, basis) == WGrad.wgrad(x^2, rshape, Mesh.interior_face, wgrad_solver)

@test wgrad_side_mon(mon_num(1), rshape, right_face, basis) == WGrad.wgrad(one, rshape, right_face, wgrad_solver)
@test wgrad_side_mon(mon_num(2), rshape, right_face, basis) == WGrad.wgrad(y, rshape, right_face, wgrad_solver)
@test_fails wgrad_side_mon(mon_num(3), rshape, right_face, basis)


# test L2 inner products of basis elements

# interior supported element inner products
int_ips = ips_interior_mons(fe_num(1), basis)
@test int_ips[1,1] == Mesh.integral_face_rel_x_face_rel_on_face(one, one, rshape, Mesh.interior_face, basis.mesh) == 1
@test int_ips[1,2] == Mesh.integral_face_rel_x_face_rel_on_face(one, y, rshape, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,3] == Mesh.integral_face_rel_x_face_rel_on_face(one, y^2, rshape, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[1,4] == Mesh.integral_face_rel_x_face_rel_on_face(one, x, rshape, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[1,5] == Mesh.integral_face_rel_x_face_rel_on_face(one, x*y, rshape, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[1,6] == Mesh.integral_face_rel_x_face_rel_on_face(one, x^2, rshape, Mesh.interior_face, basis.mesh) == 1/3

@test int_ips[2,1] == Mesh.integral_face_rel_x_face_rel_on_face(y, one, rshape, Mesh.interior_face, basis.mesh) == 1/2
@test int_ips[2,2] == Mesh.integral_face_rel_x_face_rel_on_face(y, y, rshape, Mesh.interior_face, basis.mesh) == 1/3
@test int_ips[2,3] == Mesh.integral_face_rel_x_face_rel_on_face(y, y^2, rshape, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,4] == Mesh.integral_face_rel_x_face_rel_on_face(y, x, rshape, Mesh.interior_face, basis.mesh) == 1/4
@test int_ips[2,5] == Mesh.integral_face_rel_x_face_rel_on_face(y, x*y, rshape, Mesh.interior_face, basis.mesh) == 1/6
@test int_ips[2,6] == Mesh.integral_face_rel_x_face_rel_on_face(y, x^2, rshape, Mesh.interior_face, basis.mesh) == 1/6

# side supported element inner products
top_face_ips = ips_fe_side_mons(fe_num(1), top_face, basis)

@test top_face_ips[1,1] == Mesh.integral_face_rel_x_face_rel_on_face(one, one, rshape, top_face, basis.mesh) == 1
@test top_face_ips[1,2] == Mesh.integral_face_rel_x_face_rel_on_face(one, x, rshape, top_face, basis.mesh) == 1/2
@test top_face_ips[2,1] == top_face_ips[1,2]
@test top_face_ips[2,2] == Mesh.integral_face_rel_x_face_rel_on_face(x, x, rshape, top_face, basis.mesh) == 1/3
