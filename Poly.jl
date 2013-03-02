module Poly

export Monomial,
       Polynomial,
       VectorMonomial,
       Pwr, Dim, pwr, dim, R,
       as_poly,
       domain_dim,
       monomials_of_degree_eq,
       count_monomials_of_degree_le,
       monomials_of_degree_le,
       monomial_vectors_of_degree_le,
       partial,
       divergence

import Base.string, Base.show, Base.print, Base.+, Base.*

typealias Pwr Uint8
typealias Dim Uint8
typealias R Float64

check_powers(ks...) = if any(k -> k < 0 || k > 255, ks) error("invalid power") else true end
pwr(k::Integer) = check_powers(k) && convert(Pwr,k)
dim(d::Integer) = if d < 0 || d > 255 error("invalid dimension") else convert(Dim, d) end


# A representation of a monomial function.  Information about its expression is incorporated as well
# to be used for more exact and efficient computations where possible.
type Monomial
  exps::Array{Pwr} # exponent array, length is always the dimension of the input space
  fun::Function

  Monomial(exps::Array{Pwr}) = check_powers(exps...) &&
    new(exps, mon_fun(exps...))
  Monomial(k::Integer) = check_powers(k) &&
    new([k], mon_fun(convert(Pwr,k)))
  Monomial(k1::Integer, k2::Integer) = check_powers(k1,k2) &&
    new([k1,k2], mon_fun(convert(Pwr,k1),convert(Pwr,k2)))
  Monomial(k1::Integer, k2::Integer, k3::Integer) = check_powers(k1,k2,k3) &&
    new([k1,k2,k3], mon_fun(convert(Pwr,k1),convert(Pwr,k2),convert(Pwr,k3)))
end

type Polynomial
  mons::Array{Monomial}
  coefs::Array{R}

  Polynomial(mons::Array{Monomial}, coefs::Array{R}) =
    if length(mons) != length(coefs) error("size of coefficent and monomial arrays must match")
    elseif length(mons) == 0 error("polynomial requires one or more terms")
    else
      new(mons, coefs)
    end
end

# A VectorMonomial represents a function of the form x->[0,...m(x),...0] for some monomial m,
# where all output components are zero except at m's.
type VectorMonomial
  mon::Monomial
  mon_pos::Dim # the position of the non-zero monomial component among the output components
  range_dim::Dim

  VectorMonomial(m::Monomial, mon_pos::Dim, range_dim::Dim) =
    mon_pos <= range_dim ? new(m, mon_pos, range_dim) : error("monomial component position exceeds range dimension")
end

import Base.==
==(m1::Monomial, m2::Monomial) = m1.exps == m2.exps
==(vm1::VectorMonomial, vm2::VectorMonomial) =
    vm1.mon == vm2.mon && vm1.mon_pos == vm2.mon_pos && vm1.range_dim == vm2.range_dim


# multiplication of monomials
*(m1::Monomial, m2::Monomial) = Monomial(m1.exps + m2.exps)

# addition of polynomials
+(p1::Polynomial, p2::Polynomial) = Polynomial(vcat(p1.mons, p2.mons), vcat(p1.coefs, p2.coefs))

# multiplication of polynomials
*(p1::Polynomial, p2::Polynomial) = Polynomial(
    flatten([m1 * m2 for m1 in p1.mons,  m2 in p2.mons]),
    flatten([c1 * c2 for c1 in p1.coefs, c2 in p2.coefs]))

# scalar multiplication of polynomials
*(a::R, p::Polynomial) = Polynomial(p.mons, a .* p.coefs)
*(a::Real, p::Polynomial) = Polynomial(p.mons, convert(R,a) .* p.coefs)

# mixed and promoting operations
*(c::Real, m::Monomial) = Polynomial([m],[convert(R,c)])
+(m1::Monomial, m2::Monomial) = Polynomial([m1,m2],[one(R), one(R)])
+(m::Monomial, p::Polynomial) = Polynomial(vcat(p.mons,m), vcat(p.coefs,one(R)))
+(p::Polynomial,m::Monomial) = Polynomial(vcat(p.mons,m), vcat(p.coefs,one(R)))
as_poly(m::Monomial) = Polynomial([m],[1.0])

domain_dim(m::Monomial) = convert(Dim, length(m.exps))
one_mon(dom_dim::Dim) = Monomial(zeros(Pwr, dom_dim))

mon_fun(k::Pwr) = if k == 0 x::R -> 1 else x::R -> x^k end
mon_fun(k1::Pwr, k2::Pwr) =
  if k1 == 0 && k2 == 0
    (x::R, y::R) -> 1
  elseif k1 == 0
    (x::R, y::R) -> y^k2
  elseif k2 == 0
    (x::R, y::R) -> x^k1
  else
    (x::R, y::R) -> x^k1 * y^k2
  end
mon_fun(k1::Pwr, k2::Pwr, k3::Pwr) = # (x::R, y::R, z::R) -> x^k1 * y^k2 * z^k3 # TODO: Try this version for efficiency vs the below
  if k1 == 0 && k2 == 0 && k3 == 0
    (x::R, y::R, z::R) -> 1
  elseif k1 == 0 && k2 == 0
    (x::R, y::R, z::R) -> z^k3
  elseif k1 == 0 && k3 == 0
    (x::R, y::R, z::R) -> y^k2
  elseif k2 == 0 && k3 == 0
    (x::R, y::R, z::R) -> x^k1
  elseif k1 == 0
    (x::R, y::R, z::R) -> y^k2 * z^k3
  elseif k2 == 0
    (x::R, y::R, z::R) -> x^k1 * z^k3
  elseif k3 == 0
    (x::R, y::R, z::R) -> x^k1 * y^k2
  else
    (x::R, y::R, z::R) -> x^k1 * y^k2 * z^k3
  end


# Partial derivative of a monomial.
function partial(n::Dim, m::Monomial)
  if m.exps[n] == 0
    Polynomial([one_mon(domain_dim(m))], [zero(R)])
  else
    local exps = copy(m.exps), orig_exp_n = m.exps[n]
    exps[n] = orig_exp_n - 1
    Polynomial([Monomial(exps)], [convert(R,orig_exp_n)])
  end
end

function divergence(vmon::VectorMonomial)
  partial(vmon.mon_pos, vmon.mon)
end


monomials_of_degree_eq(deg::Pwr, dom_dim::Dim) =
  if dom_dim == 1
    [Monomial(deg)]
  elseif dom_dim == 2
    [Monomial(i,deg-i) for i=zero(Pwr):deg]
  elseif dom_dim == 3
    local mons = Array(Monomial, count_monomials_of_degree_eq(deg,dom_dim)),
          pos = 1
    for i=0:deg, j=0:(deg-i)
      mons[pos] = Monomial(i,j,deg-i-j)
      pos += 1
    end
    mons
  else
    error("dimensions greater than 3 are currently not supported")
  end

count_monomials_of_degree_eq(deg::Pwr, dom_dim::Dim) = binomial(convert(Int,dom_dim + deg - 1), convert(Int,deg))
count_monomials_of_degree_le(deg::Pwr, dom_dim::Dim) = sum([count_monomials_of_degree_eq(convert(Pwr,k), dom_dim) for k=0:deg])

monomials_of_degree_le(deg::Pwr, dom_dim::Dim) =
  flatten([monomials_of_degree_eq(convert(Pwr,k), dom_dim) for k=zero(Pwr):deg]) # convert shouldn't be necessary here but currently is as of Julia 0.1.

# All MonomialVectors not exceeding the passed degree, to be used as a basis for spaces of vector valued functions with polynomial components.
monomial_vectors_of_degree_le(deg::Pwr, dom_dim::Dim, range_dim::Dim) =
  flatten(map(mon -> [VectorMonomial(mon, dim(i), dim(range_dim)) for i=1:range_dim],
              monomials_of_degree_le(deg, dom_dim)))


# String Representations

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
  string("[", "0, "^(vm.mon_pos-1), string(vm.mon), ", 0"^(vm.range_dim - vm.mon_pos), "]")
end
print(io::IO, m::VectorMonomial) = print(io, string(m))
show(io::IO, m::VectorMonomial) = print(io, m)


flatten(arrays) = vcat(arrays...)

end # end of module
