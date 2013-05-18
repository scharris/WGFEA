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
       partial,
       divergence,
       integral_on_rect_at_origin,
       drop_coefs_lt,
       coefs_closer_than,
       mons_in_same_comp_of_vmons

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

  VectorMonomial(mon::Monomial, mon_pos::Dim) =
    mon_pos <= domain_dim(mon) ? new(mon, mon_pos) :
                                 error("monomial component position exceeds monomial domain dimension")
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
  const new_coefs_by_mon = Dict{Monomial,R}()
  for i=1:length(p.mons)
    const mon = p.mons[i]
    new_coefs_by_mon[mon] = get(new_coefs_by_mon, mon, zeroR) + p.coefs[i]
  end
  mons = sort(collect(keys(new_coefs_by_mon)))
  coefs = [get(new_coefs_by_mon, mon, zeroR) for mon in mons]
  Polynomial(mons, coefs)
end

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

# real multiplication of monomials
*(a::R, m::Monomial) = Polynomial([m],[a])
*(a::Real, m::Monomial) = convert(R,a) * m

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
  sum([p.coefs[i] * reduce_dim_by_fixing(dim, val, p.mons[i]) for i in 1:length(p.coefs)])
end

function reduce_dim_by_fixing(dim::Dim, val::R, m::Monomial)
  c = val^(m.exps[dim])
  new_exps = without_dim(dim, m.exps)
  Polynomial([Monomial(new_exps)],[c])
end



# Partial derivative of a monomial.
function partial(n::Dim, m::Monomial)
  if m.exps[n] == 0
    Polynomial([one_mon(domain_dim(m))], [zeroR])
  else
    local exps = copy(m.exps), orig_exp_n = m.exps[n]
    exps[n] = orig_exp_n - 1
    Polynomial([Monomial(exps)], [convert(R,orig_exp_n)])
  end
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
  local ddim = domain_dim(m)
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
      if has(cp2_coefs_by_mon, mon)
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

one_mon(dom_dim::Dim) = Monomial(zeros(Deg, dom_dim))
zero_poly(dom_dim::Dim) = Polynomial([one_mon(dom_dim)],[zeroR])
zero_poly_vec(dom_dim::Dim) = let zerop = zero_poly(dom_dim); PolynomialVector([zerop for i=1:dom_dim]) end

flatten(arrays) = vcat(arrays...)

end # end of module
