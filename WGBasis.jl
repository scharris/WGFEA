module WGBasis

export BElNum, bel_num,
       MonNum, mon_num,
       WeakFunsPolyBasis,
       is_interior_supported,
       support_interior_num,
       is_side_supported,
       support_nb_side_num,
       interior_monomial_num,
       interior_monomial_by_num,
       interior_monomial,
       side_monomial_num,
       side_monomial_by_num,
       side_monomial,
       wgrad_for_interior_mon_num,
       wgrads_for_side_supported,
       ref_interior_mons_L2ips_matrix,
       ref_side_mons_L2ips_matrix

using Common
import Mesh, Mesh.AbstractMesh, Mesh.FENum, Mesh.NBSideInclusions, Mesh.FEFace, Mesh.fe_face
import Poly, Poly.Monomial, Poly.PolynomialVector
import WGrad, WGrad.WGradSolver

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

typealias MonNum Uint16
mon_num(i::Integer) = if i > 0 uint16(i) else error("monomial number out of range") end

type NBSideWGrads
  side_incls::NBSideInclusions
  poly_vec_on_fe1::PolynomialVector
  poly_vec_on_fe2::PolynomialVector
end

type WeakFunsPolyBasis
  interior_polys_max_deg::Deg # max degree of finite element interior polynomials
  side_polys_max_deg::Deg
  mesh::AbstractMesh

  # computed items

  dom_dim::Dim

  # the generic monomials which on interiors and sides will be used to define basis functions
  ref_interior_mons::Array{Monomial,1}
  ref_side_mons::Array{Monomial,1}
  # counts of the above
  mons_per_fe_interior::Uint
  mons_per_fe_side::Uint

  num_interior_bels::BElNum
  first_side_bel::BElNum
  total_bels::BElNum

  interior_wgrads_by_mon_num::Array{PolynomialVector,1}
  side_wgrads_by_side_mon_num::Matrix{PolynomialVector}

  # L2 inner products of basis el vs basis el for each face
  ref_interior_mons_L2ips::Matrix{R} # indexed by mon #, mon #
  ref_side_mons_L2ips::Array{R,3}    # indexed by mon #, mon #, side #

  function WeakFunsPolyBasis(interior_polys_max_deg::Deg, side_polys_max_deg::Deg, mesh::AbstractMesh)
    const dom_dim = Mesh.space_dim(mesh)
    const ref_interior_mons = Poly.monomials_of_degree_le(interior_polys_max_deg, dom_dim)
    const ref_side_mons = Poly.monomials_of_degree_le(side_polys_max_deg, dom_dim)
    const mons_per_fe_interior = length(ref_interior_mons) | uint
    const mons_per_fe_side = length(ref_side_mons) | uint
    const num_interior_bels = Mesh.num_fes(mesh) * mons_per_fe_interior

    # Precompute weak gradients by face and monomial
    const wgrad_solver = WGradSolver(deg(interior_polys_max_deg-1), mesh)

    new(interior_polys_max_deg,
        side_polys_max_deg,
        mesh,
        dom_dim,
        ref_interior_mons,
        ref_side_mons,
        mons_per_fe_interior,
        mons_per_fe_side,
        num_interior_bels,
        num_interior_bels + 1, # first_side_bel
        num_interior_bels + Mesh.num_nb_sides(mesh) * mons_per_fe_side, # total_bels
        interior_supported_bel_wgrads_by_mon_num(wgrad_solver, ref_interior_mons),
        side_supported_bel_wgrads_by_side_mon_num(wgrad_solver, Mesh.num_side_faces_per_fe(mesh), ref_side_mons),
        mon_vs_mon_L2ips_over_ref_interior(ref_interior_mons, mesh),
        mon_vs_mon_L2ips_over_ref_sides(ref_side_mons, mesh))
  end
end # ends type WeakFunsPolyBasis


function interior_supported_bel_wgrads_by_mon_num(wgrad_solver::WGradSolver, interior_mons::Array{Monomial,1})
  const wgrads = Array(PolynomialVector, length(interior_mons))
  for m=1:length(wgrads)
    wgrads[m] = WGrad.wgrad(interior_mons[m], Mesh.interior_face, wgrad_solver)
  end
  wgrads
end

function side_supported_bel_wgrads_by_side_mon_num(wgrad_solver::WGradSolver, sides_per_fe::Integer, side_mons::Array{Monomial,1})
  const mons_per_side = length(side_mons)
  const wgrads = Array(PolynomialVector, sides_per_fe, mons_per_side)
  for s=1:sides_per_fe, m=1:mons_per_side
    wgrads[s,m] = WGrad.wgrad(side_mons[m], fe_face(s), wgrad_solver)
  end
  wgrads
end

function mon_vs_mon_L2ips_over_ref_interior(mons::Array{Monomial,1}, mesh::AbstractMesh)
  const num_mons = length(mons)
  const m = Array(R, num_mons, num_mons)
  for i=1:num_mons
    const mon_i = mons[i]
    for j=1:i-1
      const ip = Mesh.integral_prod_on_ref_fe_face(mon_i, mons[j], Mesh.interior_face, mesh)
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = Mesh.integral_prod_on_ref_fe_face(mon_i, mon_i, Mesh.interior_face, mesh)
  end
  m
end

function mon_vs_mon_L2ips_over_ref_sides(mons::Array{Monomial,1}, mesh::AbstractMesh)
  const num_mons = length(mons)
  const num_sides = Mesh.num_side_faces_per_fe(mesh)
  const m = Array(R, num_mons, num_mons, num_sides)
  for s=1:num_sides
    const side = fe_face(s)
    for i=1:num_mons
      const mon_i = mons[i]
      for j=1:i-1
        const ip = Mesh.integral_prod_on_ref_fe_face(mon_i, mons[j], side, mesh)
        m[i,j,s] = ip
        m[j,i,s] = ip
      end
      m[i,i,s] = Mesh.integral_prod_on_ref_fe_face(mon_i, mon_i, side, mesh)
    end
  end
  m
end



is_interior_supported(i::BElNum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bels

function support_interior_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_interior_supported(i, basis))
  const interiors_rel_bel_ix = i-1
  div(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1
end

is_side_supported(i::BElNum, basis::WeakFunsPolyBasis) = basis.num_interior_bels < i <= basis.total_bels

function support_nb_side_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_side_supported(i, basis))
  const sides_rel_bel_ix = i - basis.first_side_bel
  div(sides_rel_bel_ix, basis.mons_per_fe_side) + 1
end

function interior_monomial_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_interior_supported(i, basis))
  const interiors_rel_bel_ix = i-1
  convert(MonNum, mod(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1)
end

function side_monomial_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_side_supported(i, basis))
  const sides_relative_bel_ix = i - basis.first_side_bel
  convert(MonNum, mod(sides_relative_bel_ix, basis.mons_per_fe_side) + 1)
end

interior_monomial_by_num(mon_num::MonNum, basis::WeakFunsPolyBasis) = basis.ref_interior_mons[mon_num]

side_monomial_by_num(mon_num::MonNum, basis::WeakFunsPolyBasis) = basis.ref_side_mons[mon_num]


interior_monomial(i::BElNum, basis::WeakFunsPolyBasis) = interior_monomial_by_num(interior_monomial_num(i, basis), basis)

side_monomial(i::BElNum, basis::WeakFunsPolyBasis) = side_monomial_by_num(side_monomial_num(i, basis), basis)

wgrad_for_interior_mon_num(mon_num::MonNum) = interior_wgrads_by_mon_num[mon_num]

function wgrads_for_side_supported(i::BElNum, basis::WeakFunsPolyBasis)
  const side_incls = Mesh.fe_inclusions_of_nb_side(i, basis)
  const mon_num = side_monomial_num(i, basis)
  const wgrad_pv_on_fe1 = side_wgrads_by_side_mon_num[side_incls.face_in_fe1, mon_num]
  const wgrad_pv_on_fe2 = side_wgrads_by_side_mon_num[side_incls.face_in_fe2, mon_num]
  NBSideWGrads(side_incls, wgrad_pv_on_fe1, wgrad_pv_on_fe2)
end

ref_interior_mons_L2ips_matrix(basis::WeakFunsPolyBasis) = basis.ref_interior_mons_L2ips

function ref_side_mons_L2ips_matrix(side::FEFace, basis::WeakFunsPolyBasis)
  if side == Mesh.interior_face error("Side face is required here.")
  else
    basis.ref_side_mons_L2ips[:,:,side]
  end
end

function isequal(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  basis1.interior_polys_max_deg == basis2.interior_polys_max_deg &&
  basis1.side_polys_max_deg == basis2.side_polys_max_deg &&
  isequal(basis1.mesh, basis2.mesh)
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
  suppt,supp,mon =
    if is_interior_supported(i, basis)
      "interior", support_interior_num(i,basis), interior_monomial_by_num(interior_monomial_num(i, basis))
    else
      begin
        incls = Mesh.fe_inclusions_of_nb_side(i, basis.mesh)
        "side", incls, side_monomial_by_num(side_monomial_num(i, basis))
      end
    end
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
    side_num = support_nb_side_num(i, basis)
    fe_incls = Mesh.fe_inclusions_of_nb_side(side_num, basis.mesh)
    [fe_incls.fe1, fe_incls.fe2]
  end
end

end # end of module
