module WGBasis

export WeakFunsPolyBasis,
       BElNum, bel_num,
       MonNum, mon_num,
       is_interior_supported,
       support_interior_num,
       is_side_supported,
       support_nb_side_num,
       fe_inclusions_of_side_support,
       ref_interior_mons,
       interior_monomial,
       interior_monomial_num,
       interior_monomial_by_num,
       ref_side_mons,
       side_monomial,
       side_monomial_num,
       side_monomial_by_face_and_num,
       mon_num_for_mon_on_side_face,
       wgrad_for_mon_on_interior,
       wgrad_for_mon_on_side_face,
       wgrads_for_side_bel,
       ips_ref_interior_mons,
       ips_ref_side_mons

using Common
import Mesh, Mesh.AbstractMesh, Mesh.FENum, Mesh.NBSideInclusions, Mesh.FERelFace, Mesh.fe_face
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
# representing the basis elements on the face are arranged lexicographically by
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
  wgrad_on_fe1::PolynomialVector
  wgrad_on_fe2::PolynomialVector
end


type WeakFunsPolyBasis

  interior_polys_max_deg::Deg # max degree of finite element interior polynomials
  side_polys_max_deg::Deg

  mesh::AbstractMesh

  # derived items

  dom_dim::Dim

  # the generic monomials which on interiors and sides will be used to define basis functions
  ref_interior_mons::Array{Monomial,1}
  ref_side_mons::Array{Array{Monomial,1},1} # indexed by side dependent dimension as declared by the mesh
  # maps to lookup monomial numbers for monomials for a given face
  # interior_mon_to_mon_num_map::Dict{Monomial, MonNum}
  mon_to_mon_num_map_by_side::Array{Dict{Monomial, MonNum},1}

  mons_per_fe_interior::Uint
  mons_per_fe_side::Uint

  # significant points in basis element enumeration
  num_interior_bels::BElNum
  first_side_bel::BElNum
  total_bels::BElNum

  # weak gradients of basis elements
  interior_mon_wgrads::Array{PolynomialVector,1} # by monomial number
  side_mon_wgrads::Matrix{PolynomialVector} # by side face, monomial number

  # L2 inner products of basis elements
  ips_ref_interior_mons::Matrix{R} # by mon #, mon #
  ips_ref_side_mons::Array{Matrix{R},1} # by side face, then mon #, mon #


  function WeakFunsPolyBasis(interior_polys_max_deg::Deg, side_polys_max_deg::Deg, mesh::AbstractMesh)
    const dom_dim = Mesh.space_dim(mesh)
    const ref_interior_mons = Poly.monomials_of_degree_le(interior_polys_max_deg, dom_dim)
    const ref_side_mons = make_ref_side_mons(side_polys_max_deg, dom_dim)
    const mons_per_fe_interior = length(ref_interior_mons) | uint
    const mons_per_fe_side = length(ref_side_mons[1]) | uint
    const num_interior_bels = Mesh.num_fes(mesh) * mons_per_fe_interior

    # Precompute weak gradients by face and monomial
    const wgrad_solver = WGradSolver(deg(interior_polys_max_deg-1), mesh)

    new(interior_polys_max_deg,
        side_polys_max_deg,
        mesh,
        dom_dim,
        ref_interior_mons,
        ref_side_mons,
        #make_interior_mon_to_mon_num_map(ref_interior_mons),
        make_side_mon_to_mon_num_maps(ref_side_mons),
        mons_per_fe_interior,
        mons_per_fe_side,
        num_interior_bels,
        num_interior_bels + 1, # first_side_bel
        num_interior_bels + Mesh.num_nb_sides(mesh) * mons_per_fe_side, # total_bels
        make_interior_supported_bel_wgrads(wgrad_solver, ref_interior_mons),
        make_side_supported_bel_wgrads(wgrad_solver, mesh, ref_side_mons),
        make_mon_vs_mon_L2_ips_over_ref_interior(ref_interior_mons, mesh),
        make_mon_vs_mon_L2_ips_over_ref_sides(ref_side_mons, mesh))
  end
end # ends type WeakFunsPolyBasis

# construction helper functions

# Side Reference Monomials
# For each dimension r, the returned array at r is the list of basis monomials with an exponent of 0
# for the r^th coordinate factor. This list will be used as a basis for those sides in the mesh within
# which dimension r is dependent (affine-dependent) on the other dimensions according to the mesh.
function make_ref_side_mons(side_polys_max_deg::Deg, dom_dim::Dim)
  const mons = Array(Array{Monomial,1}, dom_dim)
  for r=dim(1):dom_dim
    mons[r] = Poly.monomials_of_degree_le_with_const_dim(side_polys_max_deg, r, dom_dim)
  end
  mons
end

#function make_interior_mon_to_mon_num_map(ref_interior_mons::Array{Monomial,1})
#  const mon_nums_by_mon = Dict{Monomial, MonNum}()
#  const num_mons = length(ref_interior_mons)
#  sizehint(mon_nums_by_mon, num_mons)
#  for m=1:num_mons
#    mon_nums_by_mon[ref_interior_mons[m]] = m
#  end
#  mon_nums_by_mon
#end

# Return an array of monomial -> monomial # maps, as an array indexed by side dependent dimension.
function make_side_mon_to_mon_num_maps(ref_side_mons::Array{Array{Monomial,1},1})
  const d = length(ref_side_mons)
  const maps = Array(Dict{Monomial, MonNum}, d)
  for dep_dim=1:d
    const ref_mons = ref_side_mons[dep_dim]
    const num_mons = length(ref_mons)
    const mon_nums_by_mon = Dict{Monomial, MonNum}()
    sizehint(mon_nums_by_mon, num_mons)
    for m=1:num_mons
      mon_nums_by_mon[ref_mons[m]] = m
    end
    maps[dep_dim] = mon_nums_by_mon
  end
  maps
end

function make_interior_supported_bel_wgrads(wgrad_solver::WGradSolver, interior_mons::Array{Monomial,1})
  const num_mons = length(interior_mons)
  const wgrads = Array(PolynomialVector, num_mons)
  for m=1:num_mons
    wgrads[m] = WGrad.wgrad(interior_mons[m], Mesh.interior_face, wgrad_solver)
  end
  wgrads
end

function make_side_supported_bel_wgrads(wgrad_solver::WGradSolver, mesh::AbstractMesh, ref_side_mons::Array{Array{Monomial,1},1})
  const sides_per_fe = Mesh.num_side_faces_per_fe(mesh)
  const mons_per_side = length(ref_side_mons[1])
  const wgrads = Array(PolynomialVector, sides_per_fe, mons_per_side)
  for s=fe_face(1):fe_face(sides_per_fe),
      m=mon_num(1):mon_num(mons_per_side)
    const side_dep_dim = Mesh.dependent_dim_for_ref_side_face(s, mesh)
    const mon = ref_side_mons[side_dep_dim][m]
    wgrads[s,m] = WGrad.wgrad(mon, s, wgrad_solver)
  end
  wgrads
end

# Make inner products matrix for interior supported basis elements.
function make_mon_vs_mon_L2_ips_over_ref_interior(mons::Array{Monomial,1}, mesh::AbstractMesh)
  const num_mons = length(mons)
  const m = Array(R, num_mons, num_mons)
  for i=1:num_mons
    const mon_i = mons[i]
    for j=1:i-1
      const ip = Mesh.integral_face_rel_x_face_rel_on_face(mon_i, mons[j], Mesh.interior_face, mesh)
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = Mesh.integral_face_rel_x_face_rel_on_face(mon_i, mon_i, Mesh.interior_face, mesh)
  end
  m
end

# Make inner products matrices for side supported basis elements, for each side face.
function make_mon_vs_mon_L2_ips_over_ref_sides(ref_side_mons_by_dep_dim::Array{Array{Monomial,1},1}, mesh::AbstractMesh)
  const num_mons = length(ref_side_mons_by_dep_dim[1])
  const num_sides = Mesh.num_side_faces_per_fe(mesh)
  const ips = Array(Matrix{R}, num_sides)
  for s=fe_face(1):fe_face(num_sides)
    const m = Array(R, num_mons, num_mons)
    const dep_dim = Mesh.dependent_dim_for_ref_side_face(s, mesh)
    const ref_mons = ref_side_mons_by_dep_dim[dep_dim]
    for i=1:num_mons
      const mon_i = ref_mons[i]
      for j=1:i-1
        const ip = Mesh.integral_face_rel_x_face_rel_on_face(mon_i, ref_mons[j], s, mesh)
        m[i,j] = ip
        m[j,i] = ip
      end
      m[i,i] = Mesh.integral_face_rel_x_face_rel_on_face(mon_i, mon_i, s, mesh)
    end
    ips[s] = m
  end
  ips
end


# exported functions


is_interior_supported(i::BElNum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bels
is_side_supported(i::BElNum, basis::WeakFunsPolyBasis) = basis.num_interior_bels < i <= basis.total_bels

# retrieval of mesh item numbers

function support_interior_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_interior_supported(i, basis))
  const interiors_rel_bel_ix = i-1
  div(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1
end

function support_nb_side_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_side_supported(i, basis))
  const sides_rel_bel_ix = i - basis.first_side_bel
  div(sides_rel_bel_ix, basis.mons_per_fe_side) + 1
end

function fe_inclusions_of_side_support(side_beln::BElNum, basis::WeakFunsPolyBasis)
  const nb_side_num = support_nb_side_num(side_beln, basis)
  Mesh.fe_inclusions_of_nb_side(nb_side_num, basis.mesh)
end

# retrieval of monomials on support faces

ref_interior_mons(basis::WeakFunsPolyBasis) = basis.ref_interior_mons


function interior_monomial(i::BElNum, basis::WeakFunsPolyBasis)
  const interiors_rel_bel_ix = i-1
  const mon_num = mod(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1
  basis.ref_interior_mons[mon_num]
end

function interior_monomial_num(i::BElNum, basis::WeakFunsPolyBasis)
  const interiors_rel_bel_ix = i-1
  convert(MonNum, mod(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1)
end

interior_monomial_by_num(mon_num::MonNum, basis::WeakFunsPolyBasis) = basis.ref_interior_mons[mon_num]

function ref_side_mons(side_face::FERelFace, basis::WeakFunsPolyBasis)
  const side_dep_dim = Mesh.dependent_dim_for_ref_side_face(side_face, basis.mesh)
  basis.ref_side_mons[side_dep_dim]
end

function side_monomial_num(i::BElNum, basis::WeakFunsPolyBasis)
  const sides_relative_bel_ix = i - basis.first_side_bel
  convert(MonNum, mod(sides_relative_bel_ix, basis.mons_per_fe_side) + 1)
end

function side_monomial(i::BElNum, basis::WeakFunsPolyBasis)
  const nb_side_num = support_nb_side_num(i, basis)
  const side_dep_dim = Mesh.dependent_dim_for_nb_side(nb_side_num, basis.mesh)
  const mon_num = side_monomial_num(i, basis)
  basis.ref_side_mons[side_dep_dim][mon_num]
end

function side_monomial_by_face_and_num(side_face::FERelFace, mon_num::MonNum, basis::WeakFunsPolyBasis)
  const side_dep_dim = Mesh.dependent_dim_for_ref_side_face(side_face, basis.mesh)
  basis.ref_side_mons[side_dep_dim][mon_num]
end

function mon_num_for_mon_on_side_face(mon::Monomial, side_face::FERelFace, basis::WeakFunsPolyBasis)
  const side_dep_dim = Mesh.dependent_dim_for_ref_side_face(side_face, basis.mesh)
  basis.mon_to_mon_num_map_by_side[side_dep_dim][mon]
end


# weak gradient accessors

wgrad_for_mon_on_interior(mon_num::MonNum, basis::WeakFunsPolyBasis) = basis.interior_mon_wgrads[mon_num]

function wgrads_for_side_bel(i::BElNum, basis::WeakFunsPolyBasis)
  const side_incls = fe_inclusions_of_side_support(i, basis)
  const mon_num = side_monomial_num(i, basis)
  const wgrad_pv_on_fe1 = basis.side_mon_wgrads[side_incls.face_in_fe1, mon_num]
  const wgrad_pv_on_fe2 = basis.side_mon_wgrads[side_incls.face_in_fe2, mon_num]
  NBSideWGrads(side_incls, wgrad_pv_on_fe1, wgrad_pv_on_fe2)
end

function wgrad_for_mon_on_side_face(mon_num::MonNum, side_face::FERelFace, basis::WeakFunsPolyBasis)
  basis.side_mon_wgrads[side_face, mon_num]
end

# L2 inner product matrix accessors

ips_ref_interior_mons(basis::WeakFunsPolyBasis) = basis.ips_ref_interior_mons

function ips_ref_side_mons(side_face::FERelFace, basis::WeakFunsPolyBasis)
  if side_face == Mesh.interior_face error("Side face is required here, got interior.")
  else
    basis.ips_ref_side_mons[side_face]
  end
end

# ============================================
# testing and debugging aids

type BElSummary
  bel_num::BElNum
  support_type::String
  support
  mon::Monomial
end

function bel_summary(i::BElNum, basis::WeakFunsPolyBasis)
  suppt,supp,mon =
    if is_interior_supported(i, basis)
      "interior", support_interior_num(i,basis), interior_monomial(i, basis)
    else
      "side", fe_inclusions_of_side_support(i, basis), side_monomial(i, basis)
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
    const fe_incls = fe_inclusions_of_side_support(i, basis)
    [fe_incls.fe1, fe_incls.fe2]
  end
end


# equality, hashing etc

function isequal(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  basis1.interior_polys_max_deg == basis2.interior_polys_max_deg &&
  basis1.side_polys_max_deg == basis2.side_polys_max_deg &&
  isequal(basis1.mesh, basis2.mesh)
end

end # end of module
