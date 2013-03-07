module WGBasis

export BElNum, bel_num,
       WeakFunsPolyBasis,
       BElSupport,
       BElPairSupport,
       is_interior_bel,
       is_edge_bel,
       is_vert_edge_bel,
       is_horz_edge_bel,
       bel_sface_monomial,
       bel_support,
       bel_support!,
       bel_pair_support,
       bel_pair_support!,
       BElSummary,
       bel_summary,
       bel_summaries_supported_on_fe

using RMesh
using Poly


# Abbreviations & terminology
# - fe: a finite element or its number
# - finite element number: the number of the interior of the finite element
#     under the ordering described below (see "interiors ordering" diagram)
# - fes: finite elements
# - face: the interior or a single edge of a finite element
# - sface: supporting face (of a basis element), the face on which the basis element is non-zero
# - mon: monomial
# - el(s): element(s)
# - bel(s): basis element(s)
# - bel support: here, the minimal set of finite elements containing the support of the basis element.
# - interior/edge/vertical edge/horizontal edge basis element: a basis element supported
#     on a face of the type indicated (monomial on some one such face, 0 elsewhere).
# - llo: having lower left at origin. E.g. llo_rect is a rectangle having its lower left corner at the origin.
# - ix: a position number with 0 as the first position
# - rix, cix: row index, column index
# - ip: inner product, L2 unless specified otherwise.

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
# exponent, so e.g. x^0 y^2 appears prior to x^1 y^1.  Together these rules completely order
# the basis elements for our space V_h of weak function polynomials.
#
#
# With n and m being the number of columns and rows of finite elements in the mesh,
# the finite element faces are ordered as shown below:
#
# interiors ordering        vertical sides ordering       horizontal sides ordering
# ----------------------    --------------------------    ---------------------------
# |...            | mn |    |...          | m(n-1)|  |    |     |   ...   |         |
#         ...                          ...                ---------------------------
# ----------------------    --------------------------    |...            | ( m-1)n |
# | n+1 | n+2 |...| 2n |    |  n|  n+1|...| 2(n-1)|  |              ...
# ----------------------    --------------------------    ---------------------------
# |  1  |  2  |...|  n |    |  1|    2|...|  (n-1)|  |    |  1  |   ...   |    n    |
# ----------------------    --------------------------    |     |         |         |
# (also fe numbering)                                     ---------------------------


# basis element number type
typealias BElNum Uint64
bel_num(i::Integer) = if i > 0 uint64(i) else error("basis number out of range") end
const no_bel = uint64(0)

# BElSupport -  finite element faces supporting a single basis element
# This structure holds the one or two finite elements which include the support
# for a basis element, together with the number of the face which includes the
# support within each of the finite elements.  The fe1 field will always hold a
# finite element number, while for interior supported basis elements the fe2 field
# will hold RMesh.no_fe.  For vertical edges, the fe1 will contain the left finite
# element, and fe2 the right.  For horizontal edges, fe1 will contain the bottom
# finite element, and fe2 the top.  Thus the fe1 <= fe2 when fe2 is defined.
type BElSupport
  fe1::FENum
  bel_sface_fe1::FaceNum
  fe2::FENum
  bel_sface_fe2::FaceNum
end
BElSupport() = BElSupport(RMesh.no_fe, RMesh.no_face, RMesh.no_fe, RMesh.no_face)

import Base.==
function ==(bel_supp1::BElSupport, bel_supp2::BElSupport)
  bel_supp1.fe1 == bel_supp2.fe1 && bel_supp1.fe2 == bel_supp2.fe2 &&
  bel_supp1.bel_sface_fe1 == bel_supp2.bel_sface_fe1 &&
  bel_supp1.bel_sface_fe2 == bel_supp2.bel_sface_fe2
end


## BElPairSupport
## Represents the supports of a pair of basis elements, by representing
## one of three possible cases:
##   1) no finite element includes both supports:
##       => fe and face fields will be no_fe, no_face
##       => identical_multi_fe_supports will be false,
##   2) exactly one finite element includes both supports
##       => fe and face fields will hold valid fe and face numbers
##       => identical_multi_fe_supports will be false
##   3) the basis elements are supported on the same edge, so their
##      finite element/face supports (BElSupport) are identical and
##      include multiple finite elements. In this case reference the
##      BElSupport structure for the common support information.
##       => fe/face fields will be no_fe/no_face (no *unique* common fe)
##       => identical_multi_fe_supports will be true
type BElPairSupport
  fe::FENum
  bel1_sface::FaceNum
  bel2_sface::FaceNum
  identical_multi_fe_supports::Bool
end
BElPairSupport() = BElPairSupport(RMesh.no_fe, RMesh.no_face, RMesh.no_face, false)

function ==(psupp1::BElPairSupport, psupp2::BElPairSupport)
  psupp1.fe        == psupp2.fe &&
  psupp1.bel1_sface == psupp2.bel1_sface &&
  psupp1.bel2_sface == psupp2.bel2_sface
end


type WeakFunsPolyBasis
  k::Pwr # degree of finite element interior polynomials (k-1 for edges)
  dom_dim::Dim
  mesh::RectMesh

  # computed items
  fe_interior_mons::Array{Monomial,1}
  fe_edge_mons::Array{Monomial,1}
  num_fe_interior_mons::Uint
  num_fe_edge_mons::Uint
  num_interior_bels::BElNum
  num_vert_edge_bels::BElNum
  num_horz_edge_bels::BElNum
  first_vert_edge_bel::BElNum
  first_horz_edge_bel::BElNum
  total_bels::BElNum

  function WeakFunsPolyBasis(k::Pwr, dom_dim::Dim, mesh::RectMesh)
    fe_interior_mons = monomials_of_degree_le(k, dom_dim)
    fe_edge_mons = monomials_of_degree_le(pwr(k-1), dom_dim)
    num_fe_interior_mons = length(fe_interior_mons) | uint
    num_fe_edge_mons = length(fe_edge_mons) | uint
    num_interior_bels = mesh.rows * mesh.cols * num_fe_interior_mons
    num_vert_edge_bels = (mesh.cols - one(BElNum)) * mesh.rows * num_fe_edge_mons
    num_horz_edge_bels = (mesh.rows - one(BElNum)) * mesh.cols * num_fe_edge_mons
    first_vert_edge_bel = num_interior_bels + 1
    first_horz_edge_bel = num_interior_bels + num_vert_edge_bels + 1
    total_bels = num_interior_bels + num_vert_edge_bels + num_horz_edge_bels

    new(k, dom_dim, mesh,
        fe_interior_mons,
        fe_edge_mons,
        num_fe_interior_mons::Uint,
        num_fe_edge_mons::Uint,
        num_interior_bels,
        num_vert_edge_bels,
        num_horz_edge_bels,
        first_vert_edge_bel,
        first_horz_edge_bel,
        total_bels)
   end

end # ends type WeakFunsPolyBasis

function ==(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  basis1.k        == basis2.k &&
  basis1.dom_dim  == basis2.dom_dim &&
  basis1.mesh     == basis2.mesh
end


is_interior_bel(i::BElNum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bels
is_edge_bel(i::BElNum, basis::WeakFunsPolyBasis) = i > basis.num_interior_bels
is_vert_edge_bel(i::BElNum, basis::WeakFunsPolyBasis) = basis.num_interior_bels < i < basis.first_horz_edge_bel
is_horz_edge_bel(i::BElNum, basis::WeakFunsPolyBasis) = basis.first_horz_edge_bel <= i <= basis.total_bels

function bel_sface_monomial(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_bel(i, basis)
    mon_num = mod(i-1, basis.num_fe_interior_mons) + 1
    basis.fe_interior_mons[mon_num]
  elseif is_vert_edge_bel(i, basis)
    vert_edge_relative_bel_ix = i - basis.first_vert_edge_bel
    mon_num = mod(vert_edge_relative_bel_ix, basis.num_fe_edge_mons) + 1
    basis.fe_edge_mons[mon_num]
  elseif is_horz_edge_bel(i, basis)
    horz_edge_relative_bel_ix = i - basis.first_horz_edge_bel
    mon_num = mod(horz_edge_relative_bel_ix, basis.num_fe_edge_mons) + 1
    basis.fe_edge_mons[mon_num]
  else error("basis element number out of range")
  end
end

# Fill the passed BElSupport structure with information about the finite elements
# containing the support for the indicated basis element.  For vertical edges,
# the fe1 will contain the left finite element, and fe2 the right.  For horizontal
# edges, fe1 will contain the bottom finite element, and fe2 the top.
function bel_support!(i::BElNum, basis::WeakFunsPolyBasis, bel_supp::BElSupport)
  if is_interior_bel(i, basis)
    interior_num = div(i-1, basis.num_fe_interior_mons)+1 | fe_num
    bel_supp.fe1 = interior_num
    bel_supp.fe2 = RMesh.no_fe
    bel_supp.bel_sface_fe1 = RMesh.interior_face
    bel_supp.bel_sface_fe2 = RMesh.no_face
  elseif is_vert_edge_bel(i, basis)
    const mesh_cols = basis.mesh.cols
    # See "vertical sides ordering" diagram above as a reference for the vertical edge numbering.
    # First find the index (edge#-1) of the vertical edge supporting this monomial.
    vert_edge_relative_bel_ix = i - basis.first_vert_edge_bel
    edge_ix = div(vert_edge_relative_bel_ix, basis.num_fe_edge_mons)
    # Convert to edge row and column indexes (only considering the vertical edges, so rows have mesh_cols-1 items).
    edge_rix = div(edge_ix, mesh_cols-1)
    edge_cix = mod(edge_ix, mesh_cols-1)
    # Find finite element numbers on left and right of the edge (remember the fe rows have one more element than the edge rows).
    bel_supp.fe1 = edge_rix*mesh_cols + edge_cix + 1 # the fe to the left of the bel support edge
    bel_supp.fe2 = edge_rix*mesh_cols + edge_cix + 2 # the fe to the right of the bel support edge
    bel_supp.bel_sface_fe1 = RMesh.right_face # on the left fe the bel support is the right face
    bel_supp.bel_sface_fe2 = RMesh.left_face # on the right fe the bel support is the left face
  elseif is_horz_edge_bel(i, basis)
    # See "horizontal sides ordering" diagram above as a reference for the horizontal edge numbering.
    # For the horizontal edges, the edge number is the same as the number of the finite element for which it is the top edge.
    horz_edge_relative_bel_ix = i - basis.first_horz_edge_bel
    edge_num = div(horz_edge_relative_bel_ix, basis.num_fe_edge_mons) + 1
    bel_supp.fe1 = edge_num # bottom fe
    bel_supp.fe2  = edge_num + basis.mesh.cols # top fe
    bel_supp.bel_sface_fe1 = RMesh.top_face # on the bottom fe the bel support is the top face
    bel_supp.bel_sface_fe2 = RMesh.bottom_face # on the top fe the bel support is the bottom face
  else
    error("basis index out of range")
  end
end

# Functional variant of the above
bel_support(i::BElNum, basis::WeakFunsPolyBasis) =
  let bel_supp = BElSupport()
    bel_support!(i, basis, bel_supp)
    bel_supp
  end

## Returns the BElPairSupport representing the finite element supporting both
## basis elements, if any, or indicating that the basis elements are supported
## on the same edge and thus have identical multiple finite element supports.
## See the documentation for the BElPairSupport type for more details.
function bel_pair_support!(bel1_supp::BElSupport,
                           bel2_supp::BElSupport,
                           psupp::BElPairSupport)
  if bel1_supp.fe2 != RMesh.no_fe && bel1_supp == bel2_supp
    psupp.identical_multi_fe_supports = true
    psupp.fe = RMesh.no_fe
    psupp.bel1_sface = RMesh.no_face
    psupp.bel2_sface = RMesh.no_face
  elseif bel1_supp.fe1 == bel2_supp.fe1
    psupp.fe = bel1_supp.fe1
    psupp.bel1_sface = bel1_supp.bel_sface_fe1
    psupp.bel2_sface = bel2_supp.bel_sface_fe1
    psupp.identical_multi_fe_supports = false
  elseif bel1_supp.fe1 == bel2_supp.fe2
    psupp.fe = bel1_supp.fe1
    psupp.bel1_sface = bel1_supp.bel_sface_fe1
    psupp.bel2_sface = bel2_supp.bel_sface_fe2
    psupp.identical_multi_fe_supports = false
  elseif bel1_supp.fe2 == bel2_supp.fe1
    psupp.fe = bel1_supp.fe2
    psupp.bel1_sface = bel1_supp.bel_sface_fe2
    psupp.bel2_sface = bel2_supp.bel_sface_fe1
    psupp.identical_multi_fe_supports = false
  elseif bel1_supp.fe2 == bel2_supp.fe2 != RMesh.no_fe
    psupp.fe = bel1_supp.fe2
    psupp.bel1_sface = bel1_supp.bel_sface_fe2
    psupp.bel2_sface = bel2_supp.bel_sface_fe2
    psupp.identical_multi_fe_supports = false
  else # nothing in common
    psupp.fe = RMesh.no_fe
    psupp.bel1_sface = RMesh.no_face
    psupp.bel2_sface = RMesh.no_face
    psupp.identical_multi_fe_supports = false
  end
end

# Functional variant of the above
bel_pair_support(bel1_supp::BElSupport,
                 bel2_supp::BElSupport) =
  let psupp = BElPairSupport()
    bel_pair_support!(bel1_supp, bel2_supp, psupp)
    psupp
  end


# Basis Element Summary Types and Functions - should not need to use these in actual method code paths.

type BElSummary
  bel_num::BElNum
  support::BElSupport
  mon::Monomial
end

function bel_summary(i::BElNum, basis::WeakFunsPolyBasis)
  supp = bel_support(i, basis)
  mon = bel_sface_monomial(i, basis)
  BElSummary(i, supp, mon)
end


function bel_summaries_supported_on_fe(fe::FENum, basis::WeakFunsPolyBasis)
  bsums = Array(BElSummary,0)
  for i=1:basis.total_bels
    bsum = bel_summary(bel_num(i), basis)
    if bsum.support.fe1 == fe || bsum.support.fe2 == fe
        push!(bsums, bsum)
    end
  end
  bsums
end



end # end of module
