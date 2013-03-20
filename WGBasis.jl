module WGBasis

export BElNum, bel_num,
       WeakFunsPolyBasis,
       is_interior_supported,
       support_interior_num,
       is_side_supported,
       support_side_num,
       support_face_monomial

using Common
import Mesh, Mesh.AbstractMesh
import Poly, Poly.Monomial

# Abbreviations & terminology
# - fe: a finite element or its number
# - finite element number: the number of the interior of the finite element
#     under the ordering described below (see "interiors ordering" diagram)
# - fes: finite elements
# - face: the interior or a single side of a finite element
# - sface: support face (of a basis element), the face on which the basis element is non-zero
# - mon: monomial
# - el(s): element(s)
# - bel(s): basis element(s)
# - bel support: here, the minimal set of finite elements containing the support of the basis element.
# - interior/side/vertical side/horizontal side basis element: a basis element supported
#     on a face of the type indicated (monomial on some one such face, 0 elsewhere).
# - llo: having lower left at origin. E.g. llo_rect is a rectangle having its lower left corner at the origin.
# - ix: a position number with 0 as the first position
# - rix, cix: row index, column index
# - ip: inner product, L2 unless specified otherwise.
# - nb: non-boundary, not significantly intersecting the outside boundary of the mesh
#
# Let V_h^0 be the space of weak functions which are polynomials of degree k or less
# on mesh element interiors and k-1 or less on mesh element sides, and are 0 on the
# outside boundary of the mesh.
#
# Purpose: For given mesh, polynomial degree, and domain dimension, provide an object
#          representing a basis for the space V_h^0, and functions to interrogate the
#          basis to determine the supporting element, face and monomial for any given
#          basis element by number.
#
# V_h^0 Basis Elements and Numbering
#
# A basis element for V_h^0 will be a weak function which is a monomial on one face
# of a single finite element, excluding sides on the outside boundary of the mesh,
# and which is 0 on all other faces.  On interiors of the finite elements the
# maximum monomial degree will be k, and on sides k-1.
#
# The basis elements are first arranged into two groups, the first being those
# which are supported on finite element interiors, followed by those supported
# on non-boundary finite element sides. Within each of these two groups the
# faces supporting the basis elements are first ordered as determined by the
# mesh which was passed to the basis at time of construction. Finally, within
# the block of basis elements allocated to a particular face, the monomials
# representing the basis element values on the face are arranged first in
# blocks by increasing degree, then within a degree block lexicographically by
# increasing exponent, so e.g. x^0 y^2 appears prior to x^1 y^1. Together with
# the mesh which determines ordering of faces, these rules completely order the
# basis elements for our space V_h of weak function polynomials.

# basis element number type
typealias BElNum Uint64
bel_num(i::Integer) = if i > 0 uint64(i) else error("basis number out of range") end
const no_bel = uint64(0)

type WeakFunsPolyBasis
  k::Deg # max degree of finite element interior polynomials (k-1 for sides)
  dom_dim::Dim
  mesh::AbstractMesh

  # computed items

  # the generic monomials which on interiors and sides will be used to define basis functions
  per_fe_interior_mons::Array{Monomial,1}
  per_fe_side_mons::Array{Monomial,1}
  # counts of the above
  num_per_fe_interior_mons::Uint
  num_per_fe_side_mons::Uint

  num_interior_bels::BElNum
  first_side_bel::BElNum
  total_bels::BElNum

  function WeakFunsPolyBasis(k::Deg, dom_dim::Dim, mesh::AbstractMesh)
    per_fe_interior_mons = Poly.monomials_of_degree_le(k, dom_dim)
    per_fe_side_mons = Poly.monomials_of_degree_le(deg(k-1), dom_dim)
    num_per_fe_interior_mons = length(per_fe_interior_mons) | uint
    num_per_fe_side_mons = length(per_fe_side_mons) | uint
    num_interior_bels = Mesh.num_fes(mesh) * num_per_fe_interior_mons
    first_side_bel = num_interior_bels + 1
    total_bels = num_interior_bels + Mesh.num_nb_sides(mesh) * num_per_fe_side_mons
    new(k, dom_dim, mesh,
        per_fe_interior_mons,
        per_fe_side_mons,
        num_per_fe_interior_mons,
        num_per_fe_side_mons,
        num_interior_bels,
        first_side_bel,
        total_bels)
  end
end # ends type WeakFunsPolyBasis

function isequal(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  basis1.k        == basis2.k &&
  basis1.dom_dim  == basis2.dom_dim &&
  isequal(basis1.mesh, basis2.mesh)
end


is_interior_supported(i::BElNum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bels

function support_interior_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_interior_supported(i, basis))
  const interiors_rel_bel_ix = i-1
  div(interiors_rel_bel_ix, basis.num_per_fe_interior_mons) + 1
end

is_side_supported(i::BElNum, basis::WeakFunsPolyBasis) = basis.num_interior_bels < i <= basis.total_bels

function support_side_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_side_supported(i, basis))
  const sides_rel_bel_ix = i - basis.first_side_bel
  div(sides_rel_bel_ix, basis.num_per_fe_side_mons) + 1
end

function support_face_monomial(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_supported(i, basis)
    const interiors_rel_bel_ix = i-1
    const mon_num = mod(interiors_rel_bel_ix, basis.num_per_fe_interior_mons) + 1
    basis.per_fe_interior_mons[mon_num]
  else
    assert(is_side_supported(i, basis))
    const sides_relative_bel_ix = i - basis.first_side_bel
    const mon_num = mod(sides_relative_bel_ix, basis.num_per_fe_side_mons) + 1
    basis.per_fe_side_mons[mon_num]
  end
end


# ============================================
# Aids for testing and debugging

type BElSummary
  bel_num::BElNum
  support_type::String
  support
  mon::Monomial
end

function bel_summary(i::BElNum, basis::WeakFunsPolyBasis)
  suppt,supp =
    if is_interior_supported(i, basis)
      "interior", support_interior_num(i,basis)
    else
      begin
        incls = Mesh.fe_inclusions_of_nb_side(i, basis.mesh)
        "side", incls
      end
    end
  mon = support_face_monomial(i, basis)
  BElSummary(i, suppt, supp, mon)
end

function bel_summaries_supported_on_fe(fe::Mesh.FENum, basis::WeakFunsPolyBasis)
  bsums = Array(BElSummary,0)
  for i=1:basis.total_bels
    bsum = bel_summary(bel_num(i), basis)
    if bsum.support.fe1 == fe || bsum.support.fe2 == fe
        push!(bsums, bsum)
    end
  end
  bsums
end

function support_fes(i::BElNum, basis::WeakFunsPolyBasis)
  if is_interior_supported(i, basis)
    [support_interior_num(i, basis)]
  else
    side_num = support_side_num(i, basis)
    fe_incls = Mesh.fe_inclusions_of_nb_side(side_num, basis.mesh)
    [fe_incls.fe1, fe_incls.fe2]
  end
end

end # end of module
