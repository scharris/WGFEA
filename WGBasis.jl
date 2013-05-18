module WGBasis

export WeakFunsPolyBasis,
       BElNum, beln,
       MonNum, monn,
       is_interior_supported,
       ub_estimate_num_bel_bel_common_support_fe_triplets,
       support_interior_num,
       support_nb_side_num,
       fe_inclusions_of_side_support,
       mons_per_fe_interior,
       interior_mons,
       interior_mon,
       interior_mon_num,
       interior_mon_bel_num,
       fe_interior_poly_coefs,
       mons_per_fe_side,
       side_mons_for_fe_side,
       side_mons_for_oshape_side,
       side_mon_num,
       side_mon,
       side_mon_bel_num,
       fe_side_poly_coefs,
       wgrad_interior_mon,
       wgrad_side_mon,
       ips_interior_mons,
       ips_fe_side_mons,
       ips_oshape_side_mons

using Common
import Mesh, Mesh.AbstractMesh, Mesh.FENum, Mesh.NBSideInclusions, Mesh.OrientedShape, Mesh.RelFace, Mesh.rface,
       Mesh.oshape, Mesh.fen
import Poly, Poly.Monomial, Poly.Polynomial, Poly.PolynomialVector
import WGrad, WGrad.WGradSolver


# Let V_h^0 be the space of weak functions which are polynomials of degree k_i or less
# on mesh element interiors and k_s or less on mesh element sides, and are 0 on the
# outside boundary of the mesh.
#
# Purpose: For given mesh and polynomial degrees for interior and side polynomials,
#          provide an object representing a basis for the space V_h^0, and functions
#          to interrogate the basis to determine 1) the supporting finite elements,
#          relative faces and monomials for basis elements by number, 2) pre-computed
#          inner products between basis elements, and 3) pre-computed weak gradients
#          of basis elements.
#
# V_h^0 Basis Elements and Numbering
#
# A basis element for V_h^0 will be a weak function which is a monomial on one face
# of a single finite element, excluding sides on the outside boundary of the mesh,
# and which is 0 on all other faces.
#
# The basis elements are first arranged into two groups, the first being those
# which are supported on finite element interiors, followed by those supported
# on non-boundary finite element sides. The ordering of interiors and sides within
# these groups are as determined by the mesh with which the basis was constructed.
# Within the block of basis elements allocated to a particular interior or side face,
# the monomials representing the basis elements on the face are arranged
# lexicographically by increasing exponent, so e.g. x^0 y^2 appears prior to x^1 y^1.
# Together with the mesh which determines ordering of faces, these rules completely
# order the basis elements for our space V_h of weak function polynomials.

# basis element number type
typealias BElNum Uint64
beln(i::Integer) = if i > 0 convert(Uint64, i) else error("basis number out of range") end
const no_bel = uint64(0)

typealias MonNum Uint16
monn(i::Integer) = if i > 0 convert(Uint16, i) else error("monomial number out of range") end

immutable WeakFunsPolyBasis

  interior_polys_max_deg::Deg
  side_polys_max_deg::Deg

  mesh::AbstractMesh

  # derived items

  dom_dim::Dim


  # generic monomials which on interiors and sides will be used to define basis functions
  interior_mons::Array{Monomial,1}
  side_mons_by_dep_dim::Array{Array{Monomial,1},1} # indexed by side dependent dimension

  # maps to lookup monomial numbers by monomial and face
  interior_mon_to_mon_num_map::Dict{Monomial, MonNum}
  side_mon_to_mon_num_map_by_dep_dim::Array{Dict{Monomial, MonNum},1}

  mons_per_fe_interior::MonNum
  mons_per_fe_side::MonNum

  # significant points in basis element enumeration
  num_interior_bels::BElNum
  first_nb_side_bel::BElNum
  total_bels::BElNum

  wgrad_solver::WGradSolver

  # weak gradients of basis elements
  interior_mon_wgrads::Array{Array{PolynomialVector,1},1} # by fe shape, then monomial number
  side_mon_wgrads::Array{Array{Array{PolynomialVector,1},1},1} # by fe shape, then side face, then monomial number

  # L2 inner products of basis elements
  ips_interior_mons::Array{Matrix{R},1} # by fe oriented shape, then mon #, mon #
  ips_side_mons::Array{Array{Matrix{R},1},1} # by fe oriented shape, then side face, then mon #, mon #


  function WeakFunsPolyBasis(interior_polys_max_deg::Deg, side_polys_max_deg::Deg, mesh::AbstractMesh)
    const dom_dim = Mesh.space_dim(mesh)
    const interior_mons = make_interior_mons(interior_polys_max_deg, dom_dim)
    const side_mons_by_dep_dim = make_side_mons(side_polys_max_deg, dom_dim)
    const mons_per_fe_interior = monn(length(interior_mons))
    const mons_per_fe_side = monn(length(side_mons_by_dep_dim[1]))
    const num_interior_bels = Mesh.num_fes(mesh) * mons_per_fe_interior

    # Precompute weak gradients by face and monomial
    const wgrad_solver = WGradSolver(deg(interior_polys_max_deg-1), mesh)

    new(interior_polys_max_deg,
        side_polys_max_deg,
        mesh,
        dom_dim,
        interior_mons,
        side_mons_by_dep_dim,
        make_interior_mon_to_mon_num_map(interior_mons),
        make_side_mon_to_mon_num_maps(side_mons_by_dep_dim),
        mons_per_fe_interior,
        mons_per_fe_side,
        num_interior_bels,
        num_interior_bels + 1, # first_nb_side_bel
        num_interior_bels + Mesh.num_nb_sides(mesh) * mons_per_fe_side, # total_bels
        wgrad_solver,
        make_interior_mon_wgrads(interior_mons, wgrad_solver),
        make_side_mon_wgrads(side_mons_by_dep_dim, wgrad_solver),
        make_interior_mon_ips(interior_mons, mesh),
        make_side_mon_ips(side_mons_by_dep_dim, mesh))
  end
end # ends type WeakFunsPolyBasis

# equality and hashing

import Base.isequal
function isequal(basis1::WeakFunsPolyBasis, basis2::WeakFunsPolyBasis)
  if is(basis1, basis2)
    true
  else
    basis1.interior_polys_max_deg == basis2.interior_polys_max_deg &&
    basis1.side_polys_max_deg == basis2.side_polys_max_deg &&
    isequal(basis1.mesh, basis2.mesh)
  end
end

import Base.hash
hash(basis::WeakFunsPolyBasis) = basis.interior_polys_max_deg + 3 * basis.side_polys_max_deg + 5 * hash(basis.mesh)

# construction helper functions

function make_interior_mons(interior_polys_max_deg::Deg, dom_dim::Dim)
  Poly.mons_of_deg_le(interior_polys_max_deg, dom_dim)
end

# Side Reference Monomials
# For each dimension r, the returned array at r is the list of basis monomials with an exponent of 0
# for the r^th coordinate factor. This list will be used as a basis for those sides in the mesh within
# which dimension r is dependent (affine-dependent) on the other dimensions according to the mesh.
function make_side_mons(side_polys_max_deg::Deg, dom_dim::Dim)
  const mons_by_dep_dim = Array(Array{Monomial,1}, dom_dim)
  for r=dim(1):dom_dim
    mons_by_dep_dim[r] = Poly.mons_of_deg_le_with_const_dim(side_polys_max_deg, r, dom_dim)
  end
  mons_by_dep_dim
end

function make_interior_mon_to_mon_num_map(int_mons::Array{Monomial,1})
  const mon_nums_by_mon = Dict{Monomial, MonNum}()
  const num_mons = length(int_mons)
  sizehint(mon_nums_by_mon, num_mons)
  for m=1:num_mons
    mon_nums_by_mon[int_mons[m]] = m
  end
  mon_nums_by_mon
end

# Return an array of monomial -> monomial # maps, as an array indexed by side dependent dimension.
function make_side_mon_to_mon_num_maps(side_mons_by_dep_dim::Array{Array{Monomial,1},1})
  const d = length(side_mons_by_dep_dim)
  const maps = Array(Dict{Monomial, MonNum}, d)
  const num_mons = length(side_mons_by_dep_dim[1])
  for dep_dim=1:d
    const ref_mons = side_mons_by_dep_dim[dep_dim]
    const mon_nums_by_mon = Dict{Monomial, MonNum}()
    sizehint(mon_nums_by_mon, num_mons)
    for m=1:num_mons
      mon_nums_by_mon[ref_mons[m]] = m
    end
    maps[dep_dim] = mon_nums_by_mon
  end
  maps
end

# Compute weak gradients of interior supported basis element monomials, indexed by fe oriented shape,
# then monomial number.
function make_interior_mon_wgrads(int_mons::Array{Monomial,1}, wgrad_solver::WGradSolver)
  const num_oshapes = Mesh.num_oriented_element_shapes(wgrad_solver.mesh)
  const wgrads_by_oshape = Array(Array{PolynomialVector,1}, num_oshapes)
  const num_mons = length(int_mons)
  for os=oshape(1):num_oshapes
    const wgrads = Array(PolynomialVector, num_mons)
    for m=1:num_mons
      wgrads[m] = WGrad.wgrad(int_mons[m], os, Mesh.interior_face, wgrad_solver)
    end
    wgrads_by_oshape[os] = wgrads
  end
  wgrads_by_oshape
end

# Compute weak gradients of side basis elements, indexed by fe oriented shape, then by side face, then monomial number.
function make_side_mon_wgrads(side_mons_by_dep_dim::Array{Array{Monomial,1},1}, wgrad_solver::WGradSolver)
  const mesh = wgrad_solver.mesh
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  const wgrads_by_oshape = Array(Array{Array{PolynomialVector,1},1}, num_oshapes)
  const mons_per_side = length(side_mons_by_dep_dim[1])
  for os=oshape(1):num_oshapes
    const sides_per_fe = Mesh.num_side_faces_for_shape(os, mesh)
    const wgrads_by_side = Array(Array{PolynomialVector,1}, sides_per_fe)
    for sf=rface(1):sides_per_fe
      const side_dep_dim = Mesh.dependent_dim_for_oshape_side(os, sf, mesh)
      const side_mons = side_mons_by_dep_dim[side_dep_dim]
      const wgrads = Array(PolynomialVector, mons_per_side)
      for m=1:mons_per_side
        wgrads[m] = WGrad.wgrad(side_mons[m], os, sf, wgrad_solver)
      end
      wgrads_by_side[sf] = wgrads
    end
    wgrads_by_oshape[os] = wgrads_by_side
  end
  wgrads_by_oshape
end

# Make inner products matrices for interior monomials, returned in an array indexed by fe oriented shape.
function make_interior_mon_ips(mons::Array{Monomial,1}, mesh::AbstractMesh)
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  const ips_by_oshape = Array(Matrix{R}, num_oshapes)
  const num_mons = length(mons)
  for os=oshape(1):num_oshapes
    const m = Array(R, num_mons,num_mons)
    for i=1:num_mons
      const mon_i = mons[i]
      for j=1:i-1
        const ip = Mesh.integral_face_rel_x_face_rel_on_oshape_face(mon_i, mons[j], os, Mesh.interior_face, mesh)
        m[i,j] = ip
        m[j,i] = ip
      end
      m[i,i] = Mesh.integral_face_rel_x_face_rel_on_oshape_face(mon_i, mon_i, os, Mesh.interior_face, mesh)
    end
    ips_by_oshape[os] = m
  end
  ips_by_oshape
end

# Make inner products matrices for side monomials, indexed by fe oriented shape, then by side face.
function make_side_mon_ips(side_mons_by_dep_dim::Array{Array{Monomial,1},1}, mesh::AbstractMesh)
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  const ips_by_oshape = Array(Array{Matrix{R}}, num_oshapes)
  const num_mons = length(side_mons_by_dep_dim[1])
  for os=oshape(1):num_oshapes
    const num_side_faces = Mesh.num_side_faces_for_shape(os, mesh)
    const ips_by_side = Array(Matrix{R}, num_side_faces)
    for sf=rface(1):num_side_faces
      const m = Array(R, num_mons, num_mons)
      const dep_dim = Mesh.dependent_dim_for_oshape_side(os, sf, mesh)
      const ref_mons = side_mons_by_dep_dim[dep_dim]
      for i=1:num_mons
        const mon_i = ref_mons[i]
        for j=1:i-1
          const ip = Mesh.integral_face_rel_x_face_rel_on_oshape_face(mon_i, ref_mons[j], os, sf, mesh)
          m[i,j] = ip
          m[j,i] = ip
        end
        m[i,i] = Mesh.integral_face_rel_x_face_rel_on_oshape_face(mon_i, mon_i, os, sf, mesh)
      end
      ips_by_side[sf] = m
    end
    ips_by_oshape[os] = ips_by_side
  end
  ips_by_oshape
end


is_interior_supported(i::BElNum, basis::WeakFunsPolyBasis) = i <= basis.num_interior_bels

is_side_supported(i::BElNum, basis::WeakFunsPolyBasis) = basis.num_interior_bels < i <= basis.total_bels


# Returns an estimate and upper bound of the number of ordered triplets (bel1, bel2, fe) which exist
# where bel1 and bel2 are basis elements which are both supported on finite element fe.
# This function is intended to help callers allocate storage for data structures involving interacting
# basis element pairs, such as are used in the construction of sparse matrices.
function ub_estimate_num_bel_bel_common_support_fe_triplets(basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const mons_per_int = basis.mons_per_fe_interior
  const mons_per_side = basis.mons_per_fe_side
  const num_fes = Mesh.num_fes(mesh)

  sum = num_fes * mons_per_int * mons_per_int # interior mons with interior mons
  for fe=fen(1):num_fes
    const nb_sides = Mesh.num_non_boundary_sides_for_fe(fe, mesh)
    const side_sidemon_choices = nb_sides * mons_per_side
    sum += 2 * mons_per_int * nb_sides * mons_per_side + # interior mons with side mons and vice versa
           side_sidemon_choices * side_sidemon_choices  # non-boundary side mons with non-boundary side mons
  end
  sum
end

# retrieval of mesh item numbers

function support_interior_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_interior_supported(i, basis))
  const interiors_rel_bel_ix = i-1
  div(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1
end

function support_nb_side_num(i::BElNum, basis::WeakFunsPolyBasis)
  assert(is_side_supported(i, basis))
  const sides_rel_bel_ix = i - basis.first_nb_side_bel
  div(sides_rel_bel_ix, basis.mons_per_fe_side) + 1
end

function fe_inclusions_of_side_support(side_beln::BElNum, basis::WeakFunsPolyBasis)
  const nb_side_num = support_nb_side_num(side_beln, basis)
  Mesh.fe_inclusions_of_nb_side(nb_side_num, basis.mesh)
end

# retrieval of monomials on support faces

interior_mons(basis::WeakFunsPolyBasis) = basis.interior_mons

function interior_mon(i::BElNum, basis::WeakFunsPolyBasis)
  const interiors_rel_bel_ix = i-1
  const mon_num = mod(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1
  basis.interior_mons[mon_num]
end

function interior_mon_num(i::BElNum, basis::WeakFunsPolyBasis)
  const interiors_rel_bel_ix = i-1
  convert(MonNum, mod(interiors_rel_bel_ix, basis.mons_per_fe_interior) + 1)
end

function interior_mon_num(mon::Monomial, basis::WeakFunsPolyBasis)
  basis.interior_mon_to_mon_num_map[mon]
end

# TODO: unit tests
function interior_mon_bel_num(fe::FENum, monn::MonNum, basis::WeakFunsPolyBasis)
  (fe-1)*basis.mons_per_fe_interior + monn
end

# Return interior polynomial coefficients for a particular finite element given coefficients for all basis elements.
function fe_interior_poly_coefs(fe::FENum, all_basis_coefs::Vector{R}, basis::WeakFunsPolyBasis)
  const num_int_mons = basis.mons_per_fe_interior
  const first_beln = interior_mon_bel_num(fe, monn(1), basis)
  all_basis_coefs[first_beln : first_beln + num_int_mons - 1]
end

mons_per_fe_interior(basis::WeakFunsPolyBasis) = basis.mons_per_fe_interior

mons_per_fe_side(basis::WeakFunsPolyBasis) = basis.mons_per_fe_side

function side_mons_for_fe_side(fe::FENum, side_face::RelFace, basis::WeakFunsPolyBasis)
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
  const side_dep_dim = Mesh.dependent_dim_for_oshape_side(fe_oshape, side_face, basis.mesh)
  basis.side_mons_by_dep_dim[side_dep_dim]
end

function side_mons_for_oshape_side(fe_oshape::OrientedShape, side_face::RelFace, basis::WeakFunsPolyBasis)
  const side_dep_dim = Mesh.dependent_dim_for_oshape_side(fe_oshape, side_face, basis.mesh)
  basis.side_mons_by_dep_dim[side_dep_dim]
end

function side_mon_num(i::BElNum, basis::WeakFunsPolyBasis)
  const sides_relative_bel_ix = i - basis.first_nb_side_bel
  convert(MonNum, mod(sides_relative_bel_ix, basis.mons_per_fe_side) + 1)
end

function side_mon_num(mon::Monomial, fe_oshape::OrientedShape, side_face::RelFace, basis::WeakFunsPolyBasis)
  const side_dep_dim = Mesh.dependent_dim_for_oshape_side(fe_oshape, side_face, basis.mesh)
  basis.side_mon_to_mon_num_map_by_dep_dim[side_dep_dim][mon]
end


function side_mon(i::BElNum, basis::WeakFunsPolyBasis)
  const nb_side_num = support_nb_side_num(i, basis)
  const side_dep_dim = Mesh.dependent_dim_for_nb_side(nb_side_num, basis.mesh)
  const mon_num = side_mon_num(i, basis)
  basis.side_mons_by_dep_dim[side_dep_dim][mon_num]
end

function side_mon_bel_num(fe::FENum, side_face::RelFace, monn::MonNum, basis::WeakFunsPolyBasis)
  const nb_side_num = Mesh.nb_side_num_for_fe_side(fe, side_face, basis.mesh)
  basis.first_nb_side_bel + (nb_side_num-1) * basis.mons_per_fe_side + (monn-1)
end

# Return side polynomial coefficients for a particular finite element side given coefficients for all basis elements.
function fe_side_poly_coefs(fe::FENum, side_face::RelFace, all_basis_coefs::Vector{R}, basis::WeakFunsPolyBasis)
  const side_mons = side_mons_for_fe_side(fe, side_face, basis)
  const first_beln = side_mon_bel_num(fe, side_face, monn(1), basis)
  all_basis_coefs[first_beln : first_beln + length(side_mons) - 1]
end

# weak gradient accessors

wgrad_interior_mon(mon_num::MonNum, fe_oshape::OrientedShape, basis::WeakFunsPolyBasis) =
  basis.interior_mon_wgrads[fe_oshape][mon_num]

wgrad_side_mon(mon_num::MonNum, fe_oshape::OrientedShape, side_face::RelFace, basis::WeakFunsPolyBasis) =
  basis.side_mon_wgrads[fe_oshape][side_face][mon_num]

# L2 inner product matrix accessors

function ips_interior_mons(fe::FENum, basis::WeakFunsPolyBasis)
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
  basis.ips_interior_mons[fe_oshape]
end

function ips_fe_side_mons(fe::FENum, side_face::RelFace, basis::WeakFunsPolyBasis)
  if side_face == Mesh.interior_face error("Side face is required here, got interior.")
  else
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
    basis.ips_side_mons[fe_oshape][side_face]
  end
end

function ips_oshape_side_mons(fe_oshape::OrientedShape, side_face::RelFace, basis::WeakFunsPolyBasis)
  if side_face == Mesh.interior_face error("Side face is required here, got interior.")
  else
    basis.ips_side_mons[fe_oshape][side_face]
  end
end

# ============================================
# testing and debugging aids

immutable BElSummary
  beln::BElNum
  support_type::String
  support
  mon::Monomial
end

function bel_summary(i::BElNum, basis::WeakFunsPolyBasis)
  suppt,supp,mon =
    if is_interior_supported(i, basis)
      "interior", support_interior_num(i,basis), interior_mon(i, basis)
    else
      "side", fe_inclusions_of_side_support(i, basis), side_mon(i, basis)
    end
  BElSummary(i, suppt, supp, mon)
end

function bel_summaries_supported_on_fe(fe::FENum, basis::WeakFunsPolyBasis)
  bsums = Array(BElSummary,0)
  for i=1:basis.total_bels
    bsum = bel_summary(beln(i), basis)
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



end # end of module
