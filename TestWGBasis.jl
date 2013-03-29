using Test
using WGBasis

using Common
import Mesh, Mesh.fe_num
import RMesh2
import Poly, Poly.Monomial

# 2 x 3 mesh, k = 2
#  ----------
#  |  |  |  |
#  ----------
#  |  |  |  |
#  ----------
#
# 6 interior monomials, 3 side monomials
# interior basis els: 3 * 2 * 6 = 36
# vertical side basis els: 2 * 2 * 3 = 12
# horizontal side basis els: 3 * 3 = 9

basis = WeakFunsPolyBasis(deg(2), deg(1), RMesh2.RectMesh2((0.,0.), (3.,2.), 2, 3))

@test map(m -> m.exps, basis.ref_interior_mons) == {[0x00, 0x00], [0x00, 0x01], [0x01, 0x00], [0x00, 0x02], [0x01, 0x01], [0x02, 0x00]}
@test map(m -> m.exps, basis.ref_side_mons) == {[0x00, 0x00], [0x00, 0x01], [0x01, 0x00]}
@test length(basis.ref_interior_mons) == 6 == Poly.count_monomials_of_degree_le(deg(2), dim(2))
@test length(basis.ref_side_mons) == 3 == Poly.count_monomials_of_degree_le(deg(1), dim(2))
@test length(basis.ref_interior_mons) == basis.mons_per_fe_interior
@test length(basis.ref_side_mons) == basis.mons_per_fe_side

@test is_interior_supported(bel_num(1), basis)
@test is_interior_supported(bel_num(36), basis)
@test RMesh2.is_vert_nb_side(support_nb_side_num(bel_num(37), basis), basis.mesh)
@test RMesh2.is_vert_nb_side(support_nb_side_num(bel_num(48), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(49), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(55), basis), basis.mesh)
@test RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(57), basis), basis.mesh)
@test_fails !RMesh2.is_horz_nb_side(support_nb_side_num(bel_num(58), basis), basis.mesh)


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

# The next section of basis elements represent monomials on vertical sides, in blocks of 3.
@test support_fes(bel_num(37), basis) == [1,2]
@test support_fes(bel_num(38), basis) == [1,2]
@test support_fes(bel_num(39), basis) == [1,2]
@test support_fes(bel_num(40), basis) == [2,3]
@test support_fes(bel_num(41), basis) == [2,3]
@test support_fes(bel_num(42), basis) == [2,3]
@test support_fes(bel_num(43), basis) == [4,5]
@test support_fes(bel_num(44), basis) == [4,5]
@test support_fes(bel_num(45), basis) == [4,5]
@test support_fes(bel_num(46), basis) == [5,6]
@test support_fes(bel_num(47), basis) == [5,6]
@test support_fes(bel_num(48), basis) == [5,6]

# The final section of basis elements are the horizontal side monomials, assigned to sides in blocks of 3.
@test support_fes(bel_num(49), basis) == [1,4]
@test support_fes(bel_num(50), basis) == [1,4]
@test support_fes(bel_num(51), basis) == [1,4]
@test support_fes(bel_num(52), basis) == [2,5]
@test support_fes(bel_num(53), basis) == [2,5]
@test support_fes(bel_num(54), basis) == [2,5]
@test support_fes(bel_num(55), basis) == [3,6]
@test support_fes(bel_num(56), basis) == [3,6]
@test support_fes(bel_num(57), basis) == [3,6]

# Test side inclusions
incls = Mesh.NBSideInclusions()

# fe vertical side supported basis element
Mesh.fe_inclusions_of_nb_side!(support_nb_side_num(bel_num(37), basis), basis.mesh, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == RMesh2.left_face

# fe horizontal side supported basis element
Mesh.fe_inclusions_of_nb_side!(support_nb_side_num(bel_num(49), basis), basis.mesh, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == RMesh2.bottom_face

# Test support face monomial retrieval

# all interior monomials in finite element 1
@test interior_monomial(bel_num(1), basis) == Monomial(0,0)
@test interior_monomial(bel_num(2), basis) == Monomial(0,1)
@test interior_monomial(bel_num(3), basis) == Monomial(1,0)
@test interior_monomial(bel_num(4), basis) == Monomial(0,2)
@test interior_monomial(bel_num(5), basis) == Monomial(1,1)
@test interior_monomial(bel_num(6), basis) == Monomial(2,0)

# interior monomials, finite element 2
@test interior_monomial(bel_num(7), basis) == Monomial(0,0)
@test interior_monomial(bel_num(12), basis) == Monomial(2,0)

# interior monomials, finite element 6
@test interior_monomial(bel_num(31), basis) == Monomial(0,0)
@test interior_monomial(bel_num(36), basis) == Monomial(2,0)

# monomials on vertical side between finite elements 1 and 2
@test side_monomial(bel_num(37), basis) == Monomial(0,0)
@test side_monomial(bel_num(38), basis) == Monomial(0,1)
@test side_monomial(bel_num(39), basis) == Monomial(1,0)

# monomials on vertical side between finite elements 2 and 3
@test side_monomial(bel_num(40), basis) == Monomial(0,0)
@test side_monomial(bel_num(42), basis) == Monomial(1,0)

# monomials on vertical side between finite elements 5 and 6
@test side_monomial(bel_num(46), basis) == Monomial(0,0)
@test side_monomial(bel_num(47), basis) == Monomial(0,1)
@test side_monomial(bel_num(48), basis) == Monomial(1,0)

# monomials on horizontal side between finite elements 1 and 4
@test side_monomial(bel_num(49), basis) == Monomial(0,0)
@test side_monomial(bel_num(50), basis) == Monomial(0,1)
@test side_monomial(bel_num(51), basis) == Monomial(1,0)

# monomials on horizontal side between finite elements 1 and 4
@test side_monomial(bel_num(52), basis) == Monomial(0,0)
@test side_monomial(bel_num(53), basis) == Monomial(0,1)
@test side_monomial(bel_num(54), basis) == Monomial(1,0)

# monomials on horizontal side between finite elements 3 and 6
@test side_monomial(bel_num(55), basis) == Monomial(0,0)
@test side_monomial(bel_num(56), basis) == Monomial(0,1)
@test side_monomial(bel_num(57), basis) == Monomial(1,0)
