module Poly

export Monomial,
       Polynomial,
       VectorMonomial,
       PolynomialVector,
       Nomial,
       NomialOrConst,
       MonOrConst,
       canonical_form,
       value_at,
       linear_comb,
       as_poly,
       zero_poly,
       zero_poly_vec,
       domain_dim,
       count_mons_of_deg_le,
       mons_of_deg_le,
       mons_of_deg_le_with_const_dim,
       vector_mons_of_deg_le,
       reduce_dim_by_fixing,
       reduce_dim_by_subst,
       precompose_with_poly_path,
       partial,
       grad,
       antideriv,
       divergence,
       integral_on_rect_at_origin,
       translate,
       drop_coefs_lt,
       coefs_closer_than,
       mons_in_same_comp_of_vmons

require("Common.jl")
using Common


# A representation of a monomial function.  Information about its expression is incorporated as well
# to be used for more exact and efficient computations where possible.
immutable Monomial
  exps::Array{Deg,1} # exponent array, length is always the dimension of the domain

  Monomial(exps::Array{Deg,1}) = check_degs(exps...) && new(exps)
  Monomial(k::Integer) = check_degs(k) && new([convert(Deg,k)])
  Monomial(k1::Integer, k2::Integer) = check_degs(k1,k2) && new([convert(Deg,k1), convert(Deg,k2)])
  Monomial(k1::Integer, k2::Integer, k3::Integer) = check_degs(k1,k2,k3) &&
    new([convert(Deg,k1), convert(Deg,k2), convert(Deg,k3)])
  Monomial(k1::Integer, k2::Integer, k3::Integer, k4::Integer) = check_degs(k1,k2,k3,k4) &&
    new([convert(Deg,k1), convert(Deg,k2), convert(Deg,k3), convert(Deg,k4)])
  Monomial(k1::Integer, k2::Integer, k3::Integer, k4::Integer, k5::Integer) = check_degs(k1,k2,k3,k4,k5) &&
    new([convert(Deg,k1), convert(Deg,k2), convert(Deg,k3), convert(Deg,k4), convert(Deg,k5)])
end

immutable Polynomial
  mons::Array{Monomial,1}
  coefs::Array{R,1}

  Polynomial(mons::Array{Monomial,1}, coefs::Array{R,1}) =
    if length(mons) != length(coefs) error("size of coefficent and monomial arrays must match")
    elseif length(mons) == 0 error("polynomial requires one or more terms")
    else
      new(mons, coefs)
    end
end

# A VectorMonomial represents a function of the form x->[0,...m(x),...0] for some monomial m,
# where all output components are zero except at m's.
immutable VectorMonomial
  mon::Monomial
  mon_pos::Dim # the position of the non-zero monomial component among the output components

  # TODO: The new() here is causing a compiler problem, should reenable this eventually for its error checking.
  # VectorMonomial(mon::Monomial, mon_pos::Dim) =
  #   mon_pos <= domain_dim(mon) ? new(mon, mon_pos) :
  #                                error("monomial component position exceeds monomial domain dimension")
end

# A vector of polynomials having the same domain.  Or alternately it may represent a vector valued function
# where each component function is a polynomial.
immutable PolynomialVector
  polys::Array{Polynomial,1}

  function PolynomialVector(ps::Array{Polynomial,1})
    assert(length(ps) > 0, "at least one polynomial is required")
    new(ps)
  end
end

Nomial = Union(Monomial, Polynomial)
NomialOrConst = Union(Polynomial, Monomial, R)
MonOrConst = Union(Monomial, R)


import Base.isequal
isequal(m1::Monomial, m2::Monomial) = m1.exps == m2.exps
isequal(vm1::VectorMonomial, vm2::VectorMonomial) =
    vm1.mon == vm2.mon && vm1.mon_pos == vm2.mon_pos
isequal(p1::Polynomial, p2::Polynomial) =
  let cp1 = canonical_form(p1),
      cp2 = canonical_form(p2)
    cp1.mons == cp2.mons && cp1.coefs == cp2.coefs
  end
isequal(pv1::PolynomialVector, pv2::PolynomialVector) = pv1.polys == pv2.polys

import Base.hash
hash(m::Monomial) = hash(m.exps)
hash(p::Polynomial) =
  let cp = canonical_form(p)
    hash(cp.mons) + 3*hash(cp.coefs)
  end
hash(vm::VectorMonomial) = vm.mon_pos + 3*hash(vm.mon)
hash(pv::PolynomialVector) = hash(pv.polys)

import Base.isless
isless(m1::Monomial, m2::Monomial) = isless(m1.exps, m2.exps)

function canonical_form(p::Polynomial)
  # sort the monomials by exponents lexicographically, collecting like terms.
  const new_coefs_by_mon = sizehint(Dict{Monomial,R}(), length(p.mons))
  for i=1:length(p.mons)
    const mon = p.mons[i]
    new_coefs_by_mon[mon] = get(new_coefs_by_mon, mon, zeroR) + p.coefs[i]
  end

  const mons = filter(mon -> get(new_coefs_by_mon, mon, zeroR) != zeroR, sort(collect(keys(new_coefs_by_mon))))
  const coefs = [get(new_coefs_by_mon, mon, zeroR) for mon in mons]
  if length(mons) == 0
    zero_poly(domain_dim(p))
  else
    Polynomial(mons, coefs)
  end
end

# This (unfinished) method could potentially be faster (with fewer allocations) than the above.
# function canonical_form2(p::Polynomial)
#   const d = dom_dim(p)
#   const max_degs_plus1_by_dim = ones(Deg, d)
#   for mon in p.mons
#     for i=1:d
#       if mon.exps[i] + 1 > max_degs_plus1_by_dim[i]
#         max_degs_plus1_by_dim[i] = mon.exps[i] + 1
#       end
#     end
#   end

#   nzs = 0
#   const coefs_by_mon_coords = zeros(R, max_degs_plus1_by_dim...)
#   const mon_coords = zeros(Dim,d) # work array
#   for i=1:length(p.mons)
#     const mon, coef = p.mons[i], p.coefs[i]
#     # fill the coordinates array identifying this monomial
#     for j=1:d
#       mon_coords[j] = mon.exps[j] + 1
#     end
#     coefs_by_mon_coords[mon_coords...] += coef
#     if coef != zeroR nzs += 1 end
#   end

#   const mons = sizehint(Array(Monomial 0), nzs)
#   const coefs = sizehint(Array(R, 0), nzs)
#   # TODO: loop through coefs_by_mon_coords indexes, adding monomials and coefficients

#   Polynomial(mons, coefs)
# end

function canonical_form(pv::PolynomialVector)
  const n = length(pv.polys)
  polys = Array(Polynomial,n)
  for i=1:n
    polys[i] = canonical_form(pv.polys[i])
  end
  PolynomialVector(polys)
end

# multiplication operations
import Base.(*)

*(m1::Monomial, m2::Monomial) = Monomial(m1.exps + m2.exps)

*(p1::Polynomial, p2::Polynomial) = begin
  const n1 = length(p1.mons)
  const n2 = length(p2.mons)
  const n = n1 * n2
  const mons = Array(Monomial, n)
  const coefs = Array(R, n)
  i = 1
  for i1=1:n1, i2=1:n2
    mons[i] = p1.mons[i1] * p2.mons[i2]
    coefs[i] = p1.coefs[i1] * p2.coefs[i2]
    i += 1
  end
  Polynomial(mons, coefs)
end

*(m::Monomial, p::Polynomial) = begin
  const nmons = length(p.mons)
  const mons = Array(Monomial, nmons)
  for i=1:nmons
    mons[i] = m * p.mons[i]
  end
  Polynomial(mons, p.coefs)
end

*(m::Monomial, vm::VectorMonomial) =
  VectorMonomial(m * vm.mon, vm.mon_pos)
*(vm::VectorMonomial, m::Monomial) = m * vm

*(p::Polynomial, m::Monomial) = m * p

# real multiplication of polynomials
*(a::R, p::Polynomial) = Polynomial(p.mons, a .* p.coefs)
*(a::Real, p::Polynomial) = convert(R,a) * p
*(p::Polynomial, a::R) = a * p
*(p::Polynomial, a::Real) = convert(R,a) * p

# real multiplication of monomials
*(a::R, m::Monomial) = Polynomial([m],[a])
*(a::Real, m::Monomial) = convert(R,a) * m
*(m::Monomial, a::R) = a * m
*(m::Monomial, a::Real) = convert(R,a) * m

# real multiplication of PolynomialVectors
*(a::R, pv::PolynomialVector) = begin
  if a == oneR
    pv
  else
    const ncomps = length(pv.polys)
    const polys = Array(Polynomial, ncomps)
    for i=1:ncomps
      polys[i] = a * pv[i]
    end
    PolynomialVector(polys)
  end
end
*(a::Real, pv::PolynomialVector) = convert(R,a) * pv

# real multiplication of VectorMonomials
*(a::R, vm::VectorMonomial) = begin
  const dom_dim = domain_dim(vm)
  const polys = Array(Polynomial, dom_dim)
  const zeroP = zero_poly(dom_dim)
  for i=1:length(polys)
    if i == vm.mon_pos
      polys[i] = a * vm.mon
    else
      polys[i] = zeroP
    end
  end
  PolynomialVector(polys)
end
*(a::Real, vm::VectorMonomial) = convert(R,a) * vm


# dot product of vector monomials
import Base.dot
dot(vmon1::VectorMonomial, vmon2::VectorMonomial) =
  if vmon1.mon_pos == vmon2.mon_pos
    Polynomial([vmon1.mon * vmon2.mon],[oneR])
  else
    zero_poly(domain_dim(vmon1))
  end

# dot product of polynomial vectors
dot(pv1::PolynomialVector, pv2::PolynomialVector) =
  let sum = zero_poly(domain_dim(pv1[1]))
    for i=1:length(pv1.polys)
      sum += pv1.polys[i] * pv2.polys[i]
    end
    sum
  end

# addition operations
import Base.(+)

+(p1::Polynomial, p2::Polynomial) = begin
  const tot_nz = nnz(p1.coefs) + nnz(p2.coefs)
  const coefs = Array(R, tot_nz)
  const mons = Array(Monomial, tot_nz)
  nzi = 1
  for i=1:length(p1.coefs)
    if p1.coefs[i] != zeroR
      coefs[nzi] = p1.coefs[i]
      mons[nzi] = p1.mons[i]
      nzi += 1
    end
  end
  for i=1:length(p2.coefs)
    if p2.coefs[i] != zeroR
      coefs[nzi] = p2.coefs[i]
      mons[nzi] = p2.mons[i]
      nzi += 1
    end
  end
  if length(mons) > 0
    Polynomial(mons, coefs)
  else
    p1
  end
end
+(m1::Monomial, m2::Monomial) = Polynomial([m1,m2],[oneR, oneR])
+(m::Monomial, p::Polynomial) = Polynomial(vcat(p.mons,m), vcat(p.coefs, oneR))
+(p::Polynomial, m::Monomial) = Polynomial(vcat(p.mons,m), vcat(p.coefs, oneR))
+(c::R, m::Monomial) = Polynomial([one_mon(domain_dim(m)), m], [c, oneR])
+(c::Real, m::Monomial) = convert(R,c) + m
+(m::Monomial, c::R) = c + m
+(m::Monomial, c::Real) = convert(R,c) + m
+(c::R, p::Polynomial) = Polynomial(vcat(p.mons, one_mon(domain_dim(p))), vcat(p.coefs, c))
+(p::Polynomial, c::R) = c + p
+(c::Real, p::Polynomial) = convert(R,c) + p
+(p::Polynomial, c::Real) = convert(R,c) + p
+(pv1::PolynomialVector, pv2::PolynomialVector) = begin
  const ncomps = length(pv1.polys)
  assert(length(pv2.polys) == ncomps)
  const polys = Array(Polynomial, ncomps)
  for i=1:ncomps
    polys[i] = pv1.polys[i] + pv2.polys[i]
  end
  PolynomialVector(polys)
end

# minus
import Base.(-)
-(mon_or_poly_1::Nomial, mon_or_poly_2::Nomial) = mon_or_poly_1 + -mon_or_poly_2
-(m::Monomial) = Polynomial([m],[-oneR])
-(p::Polynomial) = (-oneR) * p
-(vm::VectorMonomial) = (-oneR) * vm
-(pv::PolynomialVector) = (-oneR) * pv
-(m::Monomial, c::Real) = m + (-c)
-(c::Real, m::Monomial) = c + (-m)
-(p::Polynomial, c::Real) = p + (-c)
-(c::Real, p::Polynomial) = c + (-p)

function poly_pwr(p::Polynomial, n::Deg)
  if n == 0
    one_poly(domain_dim(p))
  elseif n == 1
    p
  else
    p^n
  end
end

import Base.ref
ref(vm::VectorMonomial, i::Integer) = if i == vm.mon_pos vm.mon else zeroR end::MonOrConst
ref(pv::PolynomialVector, i::Integer) = pv.polys[i]


function linear_comb(coefs::Array{R}, vms::Array{VectorMonomial,1})
  const nmons = length(vms)
  assert(length(coefs) == nmons, "number of coeficients does not match number of vector monomials")
  assert(nmons > 1, "at least one vector monomial is required")
  sum = coefs[1] * vms[1]
  for i=2:nmons
    sum += coefs[i] * vms[i]
  end
  sum
end

as_poly(m::Monomial) = Polynomial([m],[oneR])
as_polyvec(vm::VectorMonomial) = begin
  const dom_dim = domain_dim(vm)
  const polys = Array(Polynomial, dom_dim)
  const zeroP = zero_poly(dom_dim)
  for i=1:length(polys)
    if i == vm.mon_pos
      polys[i] = vm.mon
    else
      polys[i] = zeroP
    end
  end
  PolynomialVector(polys)
end

domain_dim(m::Monomial) = convert(Dim, length(m.exps))
domain_dim(p::Polynomial) = domain_dim(p.mons[1])
domain_dim(vm::VectorMonomial) = domain_dim(vm.mon)

mons_in_same_comp_of_vmons(vmon1::VectorMonomial, vmon2::VectorMonomial) = vmon1.mon_pos == vmon2.mon_pos

function value_at(mon::Monomial, x::Vector{R})
  assert(length(mon.exps) == length(x))
  v = x[1]^mon.exps[1]
  for i=2:length(x)
    v *= x[i]^mon.exps[i]
  end
  v
end

function value_at(mon::Monomial, x::R)
  assert(length(mon.exps) == 1)
  x^(mon.exps[1])
end

function value_at(p::Polynomial, x::Vector{R})
  assert(length(p.mons[1].exps) == length(x))
  sum = zeroR
  for i in 1:length(p.mons)
    sum += p.coefs[i] * value_at(p.mons[i], x)
  end
  sum
end

function value_at(p::Polynomial, x::R)
  assert(length(p.mons[1].exps) == 1)
  sum = zeroR
  for i in 1:length(p.mons)
    sum += p.coefs[i] * value_at(p.mons[i], x)
  end
  sum
end

function value_at(pv::PolynomialVector, x::Vector{R})
  const num_comps = length(pv.polys)
  const v = zeros(R, num_comps)
  for i in 1:num_comps
    v[i] = value_at(pv.polys[i], x)
  end
  v
end

without_dim{T}(i::Integer, a::Array{T,1}) = vcat(a[1:i-1],a[i+1:])

function reduce_dim_by_fixing(dim::Dim, val::R, p::Polynomial)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * reduce_dim_by_fixing(dim, val, p.mons[i])
  end
  sum
end

function reduce_dim_by_fixing(dim::Dim, val::R, m::Monomial)
  const c = val^(m.exps[dim])
  const new_exps = without_dim(dim, m.exps)
  Polynomial([Monomial(new_exps)],[c])
end

# Substitute p in place of variable number n in monomial in_mon, yielding a new polynomial.
# The polynomial p should have one less dimension than the monomial in which it is substituted.
function reduce_dim_by_subst(n::Dim, p::Polynomial, in_mon::Monomial)
  const nth_exp = in_mon.exps[n]
  const p_to_nth_exp = poly_pwr(p, nth_exp)
  const rest_of_mon = reduce_dim_by_fixing(n, oneR, in_mon)
  p_to_nth_exp * rest_of_mon
end

# Substitute p in place of variable number n in polynomial in_poly, yielding a new polynomial.
# The polynomial p should have one less dimension than the polynomial in which is is substituted.
function reduce_dim_by_subst(n::Dim, p::Polynomial, in_poly::Polynomial)
  sum = in_poly.coefs[1] * reduce_dim_by_subst(n, p, in_poly.mons[1])
  for i=2:length(in_poly.coefs)
    sum += in_poly.coefs[i] * reduce_dim_by_subst(n, p, in_poly.mons[i])
  end
  sum
end


# Precompose the given monomial with the path through its domain space having the passed component
# polynomial functions.  The component polynomials should all be of one variable, and the number of
# components should equal the domain space dimension of the monomial.
function precompose_with_poly_path(mon::Monomial, path_comps::Array{Polynomial,1})
  const mon_dom_dim = domain_dim(mon)
  assert(mon_dom_dim == length(path_comps))

  prod = poly_pwr(path_comps[1], mon.exps[1])
  for i=2:length(path_comps)
    prod *= poly_pwr(path_comps[i], mon.exps[i])
  end
  prod
end


function precompose_with_poly_path(p::Polynomial, path_comps::Array{Polynomial,1})
  sum = p.coefs[1] * precompose_with_path(p.mons[1], path_comps)
  for i=2:length(p.coefs)
    sum += p.coefs[i] * precompose_with_path(p.mons[i], path_comps)
  end
  sum
end

# partial derivatives

function partial(n::Dim, m::Monomial)
  if m.exps[n] == 0
    Polynomial([one_mon(domain_dim(m))], [zeroR])
  else
    const exps = copy(m.exps)
    const orig_exp_n = m.exps[n]
    exps[n] = orig_exp_n - 1
    Polynomial([Monomial(exps)], [convert(R, orig_exp_n)])
  end
end

function partial(n::Dim, p::Polynomial)
  sum = p.coefs[1] * partial(dim(n), p.mons[1])
  for i=2:length(p.mons)
    sum += p.coefs[i] * partial(n, p.mons[i])
  end
  sum
end

function grad(p::Polynomial)
  const dom_dim = domain_dim(p)
  const partials = Array(Polynomial, dom_dim)
  for n=dim(1):dom_dim
    partials[n] = partial(n, p)
  end
  PolynomialVector(partials)
end

function grad(m::Monomial)
  grad(as_poly(m))
end


function divergence(vmon::VectorMonomial)
  partial(vmon.mon_pos, vmon.mon)
end


# integration

# Integrate a monomial on a rectangle of the given dimensions having lower left corner at the origin.
function integral_on_rect_at_origin(mon::Monomial, rect_dims::Array{R,1})
  assert(length(rect_dims) > 0, "rectangle dimensions must be supplied")
  prod = oneR
  for i = 1:length(rect_dims)
    deg_plus_1 = mon.exps[i] + 1
    prod *= (rect_dims[i]^deg_plus_1 / deg_plus_1)
  end
  prod
end

# Integrate a polynomial on a rectangle of the given dimensions having lower left corner at the origin.
function integral_on_rect_at_origin(p::Polynomial, rect_dims::Array{R,1})
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_on_rect_at_origin(p.mons[i], rect_dims)
  end
  sum
end

# Pre-compose the given monomial m with a translation function x -> x + v, yielding a polynomial p
# such that p(x) = m(x+v) for all x.  This allows passing points expressed relative to a new origin,
# while producing values which the original monomial would have yielded for the point relative to
# the original origin.
function translate(m::Monomial, v::Vector{R})
  const dom_dim = domain_dim(m)
  prod = oneR
  for i=dim(1):dom_dim
    prod *= poly_pwr(comp_mon(i, dom_dim) + v[i], m.exps[i])
  end
  prod
end


# Returns the polynomial function whose partial derivative in the nth input is the given monomial,
# and which has value 0 at 0.
function antideriv(n::Dim, m::Monomial)
  const exps = copy(m.exps)
  const new_exp_n = m.exps[n] + 1
  exps[n] = new_exp_n
  Polynomial([Monomial(exps)], [convert(R, 1./new_exp_n)])
end

function antideriv(n::Dim, p::Polynomial)
  sum = p.coefs[1] * antideriv(n, p.mons[1])
  for i=2:length(p.coefs)
    sum += p.coefs[i] * antideriv(n, p.mons[i])
  end
  sum
end

# Counting and listing of monomials at or not exceeding a certain degree

count_mons_of_deg_eq(deg::Deg, dom_dim::Dim) = binomial(convert(Int,dom_dim + deg - 1), convert(Int,deg))
count_mons_of_deg_le(deg::Deg, dom_dim::Dim) = sum([count_mons_of_deg_eq(convert(Deg,k), dom_dim) for k=0:deg])


function mons_of_deg_le(deg::Deg, dom_dim::Dim)
  const mons = Array(Monomial, count_mons_of_deg_le(deg, dom_dim))
  if dom_dim == 1
    pos = 1
    for i=0:deg
      mons[pos] = Monomial(i)
      pos += 1
    end
  elseif dom_dim == 2
    pos = 1
    for i=0:deg, j=0:(deg-i)
      mons[pos] = Monomial(i,j)
      pos += 1
    end
  elseif dom_dim == 3
    pos = 1
    for i=0:deg, j=0:(deg-i), k=0:(deg-i-j)
      mons[pos] = Monomial(i,j,k)
      pos += 1
    end
  elseif dom_dim == 4
    pos = 1
    for i=0:deg, j=0:(deg-i), k=0:(deg-i-j), l=0:(deg-i-j-k)
      mons[pos] = Monomial(i, j, k, l)
      pos += 1
    end
  elseif dom_dim == 5
    pos = 1
    for i=0:deg, j=0:(deg-i), k=0:(deg-i-j), l=0:(deg-i-j-k), m=0:(deg-i-j-k-l)
      mons[pos] = Monomial(i, j, k, l, m)
      pos += 1
    end
  else
    error("monomials of domain dimension $dom_dim are currently not supported")
  end
  mons
end


function mons_of_deg_le_with_const_dim(deg::Deg, const_dim::Dim, dom_dim::Dim)
  filter(m -> m.exps[const_dim] == 0,
         mons_of_deg_le(deg, dom_dim))
end

# all VectorMonomials not exceeding the indicated degree
vector_mons_of_deg_le(deg::Deg, dom_dim::Dim) =
  flatten(map(mon -> [VectorMonomial(mon, dim(i)) for i=1:domain_dim(mon)],
              mons_of_deg_le(deg, dom_dim)))



# String Representations
import Base.string, Base.show, Base.print

# string representation for monomials
function string(m::Monomial)
  const ddim = domain_dim(m)
  if ddim == 1
    "x^$(int(m.exps[1]))"
  elseif ddim == 2
    "x^$(int(m.exps[1])) y^$(int(m.exps[2]))"
  elseif ddim == 3
    "x^$(int(m.exps[1])) y^$(int(m.exps[2])) z^$(int(m.exps[3]))"
  else
    join(map(i->"x$i^$(int(m.exps[i]))", 1:ddim), " ")
  end
end
print(io::IO, m::Monomial) = print(io, string(m))
show(io::IO, m::Monomial) = print(io, m)


# string representation for polynomials
function string(p::Polynomial)
  join(map(i->"$(p.coefs[i]) $(p.mons[i])", 1:length(p.mons)), " + ")
end
print(io::IO, p::Polynomial) = print(io, string(p))
show(io::IO, p::Polynomial) = print(io, p)


# string representation for vector monomials
function string(vm::VectorMonomial)
  string("[", "0, "^(vm.mon_pos-1), string(vm.mon), ", 0"^(domain_dim(vm) - vm.mon_pos), "]")
end
print(io::IO, m::VectorMonomial) = print(io, string(m))
show(io::IO, m::VectorMonomial) = print(io, m)


# etc
function drop_coefs_lt(eps::R, p::Polynomial)
  const mons = Array(Monomial,0)
  const coefs = Array(R,0)
  for i=1:length(p.mons)
    if abs(p.coefs[i]) > eps
      push!(coefs, p.coefs[i])
      push!(mons, p.mons[i])
    end
  end
  if length(mons) != 0
    Polynomial(mons, coefs)
  else
    zero_poly(domain_dim(p))
  end
end

function drop_coefs_lt(eps::R, pv::PolynomialVector)
  const n = length(pv.polys)
  polys = Array(Polynomial,n)
  for i=1:n
    polys[i] = drop_coefs_lt(eps, pv.polys[i])
  end
  PolynomialVector(polys)
end

function coefs_by_mon(p::Polynomial)
  const cs_by_m = Dict{Monomial,R}()
  for i=1:length(p.mons)
    cs_by_m[p.mons[i]] = p.coefs[i]
  end
  cs_by_m
end

function coefs_closer_than(eps::R, p1::Polynomial, p2::Polynomial)
  const cp1 = canonical_form(p1)
  const cp2 = canonical_form(p2)
  const cp1_coefs_by_mon = coefs_by_mon(cp1)
  const cp2_coefs_by_mon = coefs_by_mon(cp2)
  for mon in keys(cp1_coefs_by_mon)
    c1 = get(cp1_coefs_by_mon, mon, zeroR)
    c2 = get(cp2_coefs_by_mon, mon, zeroR)
    if abs(c1 - c2) > eps
      return false
    else
      if haskey(cp2_coefs_by_mon, mon)
        delete!(cp2_coefs_by_mon, mon)
      end
    end
  end
  for mon in keys(cp2_coefs_by_mon)
    c1 = get(cp1_coefs_by_mon, mon, zeroR)
    c2 = get(cp2_coefs_by_mon, mon, zeroR)
    if abs(c1 - c2) > eps
      return false
    end
  end
  true
end

function coefs_closer_than(eps::R, pv1::PolynomialVector, pv2::PolynomialVector)
  assert(length(pv1.polys) == length(pv2.polys))
  for i=1:length(pv1.polys)
    if !coefs_closer_than(eps, pv1[i], pv2[i])
      return false
    end
  end
  true
end


function comp_mon(i::Dim, dom_dim::Dim)
  Monomial(one_exp_at(i, dom_dim))
end

function one_exp_at(i::Dim, dom_dim::Dim)
  const exps = zeros(Deg, dom_dim)
  exps[i] = deg(1)
  exps
end

const one_mon_1d = Monomial(0)
const one_mon_2d = Monomial(0,0)
const one_mon_3d = Monomial(0,0,0)
const one_mon_4d = Monomial(0,0,0,0)
one_mon(dom_dim::Dim) =
  if     dom_dim == 1 one_mon_1d
  elseif dom_dim == 2 one_mon_2d
  elseif dom_dim == 3 one_mon_3d
  elseif dom_dim == 4 one_mon_4d
  else Monomial(zeros(Deg, dom_dim))
  end

const one_poly_1d = Polynomial([one_mon_1d], [oneR])
const one_poly_2d = Polynomial([one_mon_2d], [oneR])
const one_poly_3d = Polynomial([one_mon_3d], [oneR])
const one_poly_4d = Polynomial([one_mon_4d], [oneR])
one_poly(dom_dim::Dim) =
  if     dom_dim == 1 one_poly_1d
  elseif dom_dim == 2 one_poly_2d
  elseif dom_dim == 3 one_poly_3d
  elseif dom_dim == 4 one_poly_4d
  else Polynomial([one_mon(dom_dim)], [oneR])
  end

const zero_poly_1d = Polynomial([one_mon_1d], [zeroR])
const zero_poly_2d = Polynomial([one_mon_2d], [zeroR])
const zero_poly_3d = Polynomial([one_mon_3d], [zeroR])
const zero_poly_4d = Polynomial([one_mon_4d], [zeroR])
zero_poly(dom_dim::Dim) =
  if     dom_dim == 1 zero_poly_1d
  elseif dom_dim == 2 zero_poly_2d
  elseif dom_dim == 3 zero_poly_3d
  elseif dom_dim == 4 zero_poly_4d
  else Polynomial([one_mon(dom_dim)], [zeroR])
  end

const zero_polyvec_1d = PolynomialVector([zero_poly_1d])
const zero_polyvec_2d = PolynomialVector([zero_poly_2d for i=1:2])
const zero_polyvec_3d = PolynomialVector([zero_poly_3d for i=1:3])
const zero_polyvec_4d = PolynomialVector([zero_poly_4d for i=1:4])
zero_poly_vec(dom_dim::Dim) =
  if     dom_dim == 1 zero_polyvec_1d
  elseif dom_dim == 2 zero_polyvec_2d
  elseif dom_dim == 3 zero_polyvec_3d
  elseif dom_dim == 4 zero_polyvec_4d
  else
    let zerop = zero_poly(dom_dim); PolynomialVector([zerop for i=1:dom_dim]) end
  end

flatten(arrays) = vcat(arrays...)

end # end of module
