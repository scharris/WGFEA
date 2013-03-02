module WGBasis

export BENum, be_num,
       WeakFunsPolyBasis,
       BESupport,
       BEPairSupport,
       is_interior_be,
       is_edge_be,
       is_vert_edge_be,
       is_horz_edge_be,
       be_sface_monomial,
       be_support,
       be_support!,
       be_pair_support,
       be_pair_support!

using RMesh
using Poly


# Abbreviations & terminology
# - fe: a finite element or its number
# - finite element number: the number of the interior of the finite element
#     under the ordering described below (see "interiors ordering" diagram)
# - fes: finite elements
# - face: the interior or a single edge of a finite element
# - sface: supporting face (of a basis element), the face on which the be is non-zero
# - mon: monomial
# - el(s): element(s)
# - be(s): basis element(s)
# - be support - here, the minimal set of finite elements containing the support of the basis element.
# - interior/edge/vertical edge/horizontal edge basis element: a basis element supported
#   on a face of the type indicated (monomial on some one such face, 0 elsewhere).
# - ix: a position number with 0 as the first position
# - rix, cix: row index, column index


# V_h Basis Elements and Numbering
#
# A basis element for V_h will be a weak function which is a monomial on one face
# of a single finite element, excluding edges on the outside boundary of the mesh,
# and which is 0 on all other faces.  On interiors of the finite elements the
# maximum monomial degree will be k, and on edges k-1.
#
# The basis elements are arranged first by the type of face on which they are non-zero,
# with those which are non-zero on interiors first, then those on vertical edges, and
# finally those on horizontal edges.  Within these three groups the elements are
# arranged by increasing finite element face numbers, with the faces of a particular
# type ordered starting at the lower left of the mesh and proceeding to the right on
# the same row, begining again on the next row at the left and so on (see diagrams below).
# Within a block of basis elements for a particular face, the monomials are arranged first
# in blocks by increasing degree, then within a degree block lexicographically by increasing
# exponent, so e.g. x^0 y^2 appears prior to x^1 y^1.
#
# With n and m being the number of columns and rows of finite elements in the mesh,
# the finite element faces are ordered as shown below:
#
# interiors ordering        vertical sides ordering       horizontal sides ordering
# ----------------------    --------------------------    --------------------------
# |...            | mn |    |...          | m(n-1)|  |    |     |   ...   |         |
#         ...                          ...                ---------------------------
# ----------------------    --------------------------    |...            | ( m-1)n |
# | n+1 | n+2 |...| 2n |    |  n|  n+1|...| 2(n-1)|  |              ...
# ----------------------    --------------------------    ---------------------------
# |  1  |  2  |...|  n |    |  1|    2|...|  (n-1)|  |    |  1  |   ...   |    n    |
# ----------------------    --------------------------    |     |         |         |
#                                                         ---------------------------

# basis element number type
typealias BENum Uint64
be_num(i::Integer) = if i > 0 uint64(i) else error("basis number out of range") end
const no_be = uint64(0)

# BESupport -  finite element faces supporting a single basis element
# This structure holds the one or two finite elements including the support
# for a basis element, together with the face number including the support
# in each of the finite elements.  The fe1 field will always hold a finite
# element number, while for interior supported basis elements the fe2 field
# will hold RMesh.no_fe.  For vertical edges, the fe1 will contain the left finite
# element, and fe2 the right.  For horizontal edges, fe1 will contain the bottom
# finite element, and fe2 the top.  Thus the fe1 <= fe2 when fe2 is defined.
type BESupport
  fe1::FENum
  be_sface_fe1::FaceNum
  fe2::FENum
  be_sface_fe2::FaceNum
end
BESupport() = BESupport(RMesh.no_fe, RMesh.no_face, RMesh.no_fe, RMesh.no_face)

import Base.==
function ==(be_supp1::BESupport, be_supp2::BESupport)
  be_supp1.fe1 == be_supp2.fe1 && be_supp1.fe2 == be_supp2.fe2 &&
  be_supp1.be_sface_fe1 == be_supp2.be_sface_fe1 &&
  be_supp1.be_sface_fe2 == be_supp2.be_sface_fe2
end


## BEPairSupport
## Represents the supports of a pair of basis elements, by representing
## one of three possible cases:
##   1) no finite element includes both supports:
##       => fe and face fields are no_fe, no_face
##       => identical_multi_fe_supports is false,
##   2) exactly one finite element includes both supports
##       => fe and face fields hold fe and face numbers
##       => identical_multi_fe_supports is false
##   3) the basis elements are supported on the same edge, so their
##      finite element/face supports (BESupport) are identical and
##      include multiple finite elements. In this case reference the
##      BESupport structure for the common support information.
##       => fe/face fields = no_fe/no_face (no *unique* common fe)
##       => identical_multi_fe_supports is true
type BEPairSupport
  fe::FENum
  be1_sface::FaceNum
  be2_sface::FaceNum
  identical_multi_fe_supports::Bool
end
BEPairSupport() = BEPairSupport(RMesh.no_fe, RMesh.no_face, RMesh.no_face, false)

function ==(psupp1::BEPairSupport, psupp2::BEPairSupport)
  psupp1.fe        == psupp2.fe &&
  psupp1.be1_sface == psupp2.be1_sface &&
  psupp1.be2_sface == psupp2.be2_sface
end


type WeakFunsPolyBasis
  k::Pwr # degree of finite element interior polynomials (k-1 for edges)
  dom_dim::Dim
  mesh::RectMesh

  # computed items
  fe_interior_mons::Array{Monomial}
  fe_edge_mons::Array{Monomial}
  num_fe_interior_mons::Uint
  num_fe_edge_mons::Uint
  num_interior_bes::BENum
  num_vert_edge_bes::BENum
  num_horz_edge_bes::BENum
  first_vert_edge_be::BENum
  first_horz_edge_be::BENum
  total_bes::BENum

  function WeakFunsPolyBasis(k::Pwr, dom_dim::Dim, mesh::RectMesh)
    fe_interior_mons = monomials_of_degree_le(k, dom_dim)
    fe_edge_mons = monomials_of_degree_le(pwr(k-1), dom_dim)
    num_fe_interior_mons = length(fe_interior_mons) | uint
    num_fe_edge_mons = length(fe_edge_mons) | uint
    num_interior_bes = mesh.rows * mesh.cols * num_fe_interior_mons
    num_vert_edge_bes = (mesh.cols - one(BENum)) * mesh.rows * num_fe_edge_mons
    num_horz_edge_bes = (mesh.rows - one(BENum)) * mesh.cols * num_fe_edge_mons
    first_vert_edge_be = num_interior_bes + 1
    first_horz_edge_be = num_interior_bes + num_vert_edge_bes + 1
    total_bes = num_interior_bes + num_vert_edge_bes + num_horz_edge_bes

    new(k, dom_dim, mesh,
        fe_interior_mons,
        fe_edge_mons,
        num_fe_interior_mons::Uint,
        num_fe_edge_mons::Uint,
        num_interior_bes,
        num_vert_edge_bes,
        num_horz_edge_bes,
        first_vert_edge_be,
        first_horz_edge_be,
        total_bes)
   end

end # ends type WeakFunsPolyBasis

function ==(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  basis1.k        == basis2.k &&
  basis1.dom_dim  == basis2.dom_dim &&
  basis1.mesh     == basis2.mesh
end


is_interior_be(i::BENum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bes
is_edge_be(i::BENum, basis::WeakFunsPolyBasis) = i > basis.num_interior_bes
is_vert_edge_be(i::BENum, basis::WeakFunsPolyBasis) = basis.num_interior_bes < i < basis.first_horz_edge_be
is_horz_edge_be(i::BENum, basis::WeakFunsPolyBasis) = basis.first_horz_edge_be <= i <= basis.total_bes

function be_sface_monomial(i::BENum, basis::WeakFunsPolyBasis)
  if is_interior_be(i, basis)
    mon_num = mod(i-1, basis.num_fe_interior_mons) + 1
    basis.fe_interior_mons[mon_num]
  elseif is_vert_edge_be(i, basis)
    vert_edge_relative_be_ix = i - basis.first_vert_edge_be
    mon_num = mod(vert_edge_relative_be_ix, basis.num_fe_edge_mons) + 1
    basis.fe_edge_mons[mon_num]
  elseif is_horz_edge_be(i, basis)
    horz_edge_relative_be_ix = i - basis.first_horz_edge_be
    mon_num = mod(horz_edge_relative_be_ix, basis.num_fe_edge_mons) + 1
    basis.fe_edge_mons[mon_num]
  else error("basis element number out of range")
  end
end

# Fill the passed BESupport structure with information about the finite elements
# containing the support for the indicated basis element.  For vertical edges,
# the fe1 will contain the left finite element, and fe2 the right.  For horizontal
# edges, fe1 will contain the bottom finite element, and fe2 the top.
function be_support!(i::BENum, basis::WeakFunsPolyBasis, be_supp::BESupport)
  if is_interior_be(i, basis)
    interior_num = div(i-1, basis.num_fe_interior_mons)+1 | fe_num
    be_supp.fe1 = interior_num
    be_supp.fe2 = RMesh.no_fe
    be_supp.be_sface_fe1 = RMesh.interior_face
    be_supp.be_sface_fe2 = RMesh.no_face
  elseif is_vert_edge_be(i, basis)
    const mesh_cols = basis.mesh.cols
    # See "vertical sides ordering" diagram above as a reference for the vertical edge numbering.
    # First find the index (edge#-1) of the vertical edge supporting this monomial.
    vert_edge_relative_be_ix = i - basis.first_vert_edge_be
    edge_ix = div(vert_edge_relative_be_ix, basis.num_fe_edge_mons)
    # Convert to edge row and column indexes (only considering the vertical edges, so rows have mesh_cols-1 items).
    edge_rix = div(edge_ix, mesh_cols-1)
    edge_cix = mod(edge_ix, mesh_cols-1)
    # Find finite element numbers on left and right of the edge (remember the fe rows have one more element than the edge rows).
    be_supp.fe1 = edge_rix*mesh_cols + edge_cix + 1 # the fe to the left of the be support edge
    be_supp.fe2 = edge_rix*mesh_cols + edge_cix + 2 # the fe to the right of the be support edge
    be_supp.be_sface_fe1 = RMesh.right_face # on the left fe the be support is the right face
    be_supp.be_sface_fe2 = RMesh.left_face # on the right fe the be support is the left face
  elseif is_horz_edge_be(i, basis)
    # See "horizontal sides ordering" diagram above as a reference for the horizontal edge numbering.
    # For the horizontal edges, the edge number is the same as the number of the finite element for which it is the top edge.
    horz_edge_relative_be_ix = i - basis.first_horz_edge_be
    edge_num = div(horz_edge_relative_be_ix, basis.num_fe_edge_mons) + 1
    be_supp.fe1 = edge_num # bottom fe
    be_supp.fe2  = edge_num + basis.mesh.cols # top fe
    be_supp.be_sface_fe1 = RMesh.top_face # on the bottom fe the be support is the top face
    be_supp.be_sface_fe2 = RMesh.bottom_face # on the top fe the be support is the bottom face
  else
    error("basis index out of range")
  end
end

# Functional variant of the above
be_support(i::BENum, basis::WeakFunsPolyBasis) =
  let be_supp = BESupport()
    be_support!(i, basis, be_supp)
    be_supp
  end

## Returns the BEPairSupport representing the finite element supporting both
## basis elements, if any, or indicating that the basis elements are supported
## on the same edge and thus have identical multiple finite element supports.
## See the documentation for the BEPairSupport type for more details.
function be_pair_support!(be1_supp::BESupport,
                          be2_supp::BESupport,
                          psupp::BEPairSupport)
  if be1_supp.fe2 != RMesh.no_fe && be1_supp == be2_supp
    psupp.identical_multi_fe_supports = true
    psupp.fe = RMesh.no_fe
    psupp.be1_sface = RMesh.no_face
    psupp.be2_sface = RMesh.no_face
  elseif be1_supp.fe1 == be2_supp.fe1
    psupp.fe = be1_supp.fe1
    psupp.be1_sface = be1_supp.be_sface_fe1
    psupp.be2_sface = be2_supp.be_sface_fe1
    psupp.identical_multi_fe_supports = false
  elseif be1_supp.fe1 == be2_supp.fe2
    psupp.fe = be1_supp.fe1
    psupp.be1_sface = be1_supp.be_sface_fe1
    psupp.be2_sface = be2_supp.be_sface_fe2
    psupp.identical_multi_fe_supports = false
  elseif be1_supp.fe2 == be2_supp.fe1
    psupp.fe = be1_supp.fe2
    psupp.be1_sface = be1_supp.be_sface_fe2
    psupp.be2_sface = be2_supp.be_sface_fe1
    psupp.identical_multi_fe_supports = false
  elseif be1_supp.fe2 == be2_supp.fe2 != RMesh.no_fe
    psupp.fe = be1_supp.fe2
    psupp.be1_sface = be1_supp.be_sface_fe2
    psupp.be2_sface = be2_supp.be_sface_fe2
    psupp.identical_multi_fe_supports = false
  else # nothing in common
    psupp.fe = RMesh.no_fe
    psupp.be1_sface = RMesh.no_face
    psupp.be2_sface = RMesh.no_face
    psupp.identical_multi_fe_supports = false
  end
end

# Functional variant of the above
be_pair_support(be1_supp::BESupport,
                be2_supp::BESupport) =
  let psupp = BEPairSupport()
    be_pair_support!(be1_supp, be2_supp, psupp)
    psupp
  end


end # end of module
