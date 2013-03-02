using WGBasis
using RMesh
using Poly

# 3 cols x 2 rows mesh, k = 2
#  ----------
#  |  |  |  |
#  ----------
#  |  |  |  |
#  ----------
#
# 6 interior monomials, 3 edge monomials
# interior basis els: 3 * 2 * 6 = 36
# vertical edge basis els: 2 * 2 * 3 = 12
# horizontal edge basis els: 3 * 3 = 9

rmesh = RectMesh((0.,0.), (3.,2.), fe_coord(2), fe_coord(3))
basis = WeakFunsPolyBasis(pwr(2), dim(2), rmesh)

assert(map(m -> m.exps, basis.fe_interior_mons) == {[0x00, 0x00], [0x00, 0x01], [0x01, 0x00], [0x00, 0x02], [0x01, 0x01], [0x02, 0x00]})
assert(map(m -> m.exps, basis.fe_edge_mons) == {[0x00, 0x00], [0x00, 0x01], [0x01, 0x00]})
assert(length(basis.fe_interior_mons) == 6 == count_monomials_of_degree_le(pwr(2), dim(2)))
assert(length(basis.fe_edge_mons) == 3 == count_monomials_of_degree_le(pwr(1), dim(2)))

assert(is_interior_be(fe_num(1), basis))
assert(is_interior_be(fe_num(36), basis))
assert(is_vert_edge_be(fe_num(37), basis))
assert(is_vert_edge_be(fe_num(48), basis))
assert(is_horz_edge_be(fe_num(49), basis))
assert(is_horz_edge_be(fe_num(55), basis))
assert(is_horz_edge_be(fe_num(57), basis))
assert(!is_horz_edge_be(fe_num(58), basis))


be_support_fes = (i::BENum, basis::WeakFunsPolyBasis) ->
  let be_supp = be_support(i, basis)
    if be_supp.fe2 == RMesh.no_fe
      [be_supp.fe1]
    else
      [be_supp.fe1, be_supp.fe2]
    end
  end

# The first section of basis elements are monomials which are assigned to interiors in blocks of 6.
assert(be_support_fes(be_num(1), basis) == [1])
assert(be_support_fes(be_num(6), basis) == [1])
assert(be_support_fes(be_num(7), basis) == [2])
assert(be_support_fes(be_num(12), basis) == [2])
assert(be_support_fes(be_num(13), basis) == [3])
assert(be_support_fes(be_num(30), basis) == [5])
assert(be_support_fes(be_num(31), basis) == [6])
assert(be_support_fes(be_num(36), basis) == [6])

# The next section of basis elements represent monomials on vertical edges, in blocks of 3.
assert(be_support_fes(be_num(37), basis) == [1,2])
assert(be_support_fes(be_num(38), basis) == [1,2])
assert(be_support_fes(be_num(39), basis) == [1,2])
assert(be_support_fes(be_num(40), basis) == [2,3])
assert(be_support_fes(be_num(41), basis) == [2,3])
assert(be_support_fes(be_num(42), basis) == [2,3])
assert(be_support_fes(be_num(43), basis) == [4,5])
assert(be_support_fes(be_num(44), basis) == [4,5])
assert(be_support_fes(be_num(45), basis) == [4,5])
assert(be_support_fes(be_num(46), basis) == [5,6])
assert(be_support_fes(be_num(47), basis) == [5,6])
assert(be_support_fes(be_num(48), basis) == [5,6])

# The final section of basis elements are the horizontal edge monomials, assigned to edges in blocks of 3.
assert(be_support_fes(be_num(49), basis) == [1,4])
assert(be_support_fes(be_num(50), basis) == [1,4])
assert(be_support_fes(be_num(51), basis) == [1,4])
assert(be_support_fes(be_num(52), basis) == [2,5])
assert(be_support_fes(be_num(53), basis) == [2,5])
assert(be_support_fes(be_num(54), basis) == [2,5])
assert(be_support_fes(be_num(55), basis) == [3,6])
assert(be_support_fes(be_num(56), basis) == [3,6])
assert(be_support_fes(be_num(57), basis) == [3,6])


# Test the more specific structure-filling function for finding the support of a basis element.
be_supp = BESupport()

# fe interior supported basis element
be_support!(be_num(1), basis, be_supp)
assert(be_supp.fe1 == fe_num(1))
assert(be_supp.fe2 == RMesh.no_fe)
assert(be_supp.be_sface_fe1 == RMesh.interior_face)

# fe vertical edge supported basis element
be_support!(be_num(37), basis, be_supp)
assert(be_supp.fe1 == fe_num(1))
assert(be_supp.be_sface_fe1 == RMesh.right_face)
assert(be_supp.fe2 == fe_num(2))
assert(be_supp.be_sface_fe2 == RMesh.left_face)

# fe horizontal edge supported basis element
be_support!(be_num(49), basis, be_supp)
assert(be_supp.fe1 == fe_num(1))
assert(be_supp.be_sface_fe1 == RMesh.top_face)
assert(be_supp.fe2 == fe_num(4))
assert(be_supp.be_sface_fe2 == RMesh.bottom_face)


#
# Test finding finite element containing supports of both basis elements of a pair.
#

function test_pair_support(be1::BENum, be2::BENum,
                           expected_fe::FENum, expected_face1::FaceNum, expected_face2::FaceNum, expected_identical_multi_fe_supports::Bool)
  be_supp1 = BESupport()
  be_supp2 = BESupport()
  be_support!(be1, basis, be_supp1)
  be_support!(be2, basis, be_supp2)
  pair_supp = be_pair_support(be_supp1, be_supp2)
  assert(pair_supp.fe == expected_fe)
  assert(pair_supp.be1_sface == expected_face1)
  assert(pair_supp.be2_sface == expected_face2)
  assert(pair_supp.identical_multi_fe_supports == expected_identical_multi_fe_supports)
end

# common support finite element for pair of basis elements supported on same interior
test_pair_support(be_num(1), be_num(1), fe_num(1), RMesh.interior_face, RMesh.interior_face, false)

# both basis elements are still supported on interior 1
test_pair_support(be_num(1), be_num(6), fe_num(1), RMesh.interior_face, RMesh.interior_face, false)

# second basis element now supported on a different interior, leaving no common supporting fe.
test_pair_support(be_num(1), be_num(7), RMesh.no_fe, RMesh.no_face, RMesh.no_face, false)

# second basis element now is the right edge of fe 1
test_pair_support(be_num(1), be_num(37), fe_num(1), RMesh.interior_face, RMesh.right_face, false)
test_pair_support(be_num(1), be_num(39), fe_num(1), RMesh.interior_face, RMesh.right_face, false)

# second edge is now the next vertical edge over, no longer supported in fe 1, leaving no common supporting fe.
test_pair_support(be_num(1), be_num(40),  RMesh.no_fe, RMesh.no_face, RMesh.no_face, false)

# test an interior be of fe 1 against a horizontal edge be of fe 1
test_pair_support(be_num(1), be_num(51), fe_num(1), RMesh.interior_face, RMesh.top_face, false)

# the next basis element is on the next horizontal edge over, leaving no common supporting fe.
test_pair_support(be_num(1), be_num(52), RMesh.no_fe, RMesh.no_face, RMesh.no_face, false)


# test a basis element supported on the vertical edge between fes 2 and 3, with another supported on the same edge.
test_pair_support(be_num(41), be_num(42), RMesh.no_fe, RMesh.no_face, RMesh.no_face, true)

# Test supporting face monomial retrieval

# all interior monomials in finite element 1
assert(be_sface_monomial(be_num(1), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(2), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(3), basis) == Monomial(1,0))
assert(be_sface_monomial(be_num(4), basis) == Monomial(0,2))
assert(be_sface_monomial(be_num(5), basis) == Monomial(1,1))
assert(be_sface_monomial(be_num(6), basis) == Monomial(2,0))

# interior monomials, finite element 2
assert(be_sface_monomial(be_num(7), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(12), basis) == Monomial(2,0))

# interior monomials, finite element 6
assert(be_sface_monomial(be_num(31), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(36), basis) == Monomial(2,0))

# monomials on vertical edge between finite elements 1 and 2
assert(be_sface_monomial(be_num(37), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(38), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(39), basis) == Monomial(1,0))

# monomials on vertical edge between finite elements 2 and 3
assert(be_sface_monomial(be_num(40), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(42), basis) == Monomial(1,0))

# monomials on vertical edge between finite elements 5 and 6
assert(be_sface_monomial(be_num(46), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(47), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(48), basis) == Monomial(1,0))

# monomials on horizontal edge between finite elements 1 and 4
assert(be_sface_monomial(be_num(49), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(50), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(51), basis) == Monomial(1,0))

# monomials on horizontal edge between finite elements 1 and 4
assert(be_sface_monomial(be_num(52), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(53), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(54), basis) == Monomial(1,0))

# monomials on horizontal edge between finite elements 3 and 6
assert(be_sface_monomial(be_num(55), basis) == Monomial(0,0))
assert(be_sface_monomial(be_num(56), basis) == Monomial(0,1))
assert(be_sface_monomial(be_num(57), basis) == Monomial(1,0))
