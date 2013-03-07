using Poly

x = Monomial(1,0)
y = Monomial(0,1)

assert(x * y == Monomial(1,1))
assert(x * y^2 == Monomial(1,2))
assert(domain_dim(x) == 2)

assert(monomial_value(x * y^2, 1., 2.) == 4.)
assert(polynomial_value(2x + 2(x*y), 1.0, 2.0) == 6.)


function test_poly_vals_eq(p1::Polynomial, p2::Polynomial)
  for x in [0.0, 1.0, 2.0, 3.0, 4.0], y in [10.0, 11.0, 12.0, 13.0]
    assert(polynomial_value(p1, x, y) == polynomial_value(p2, x, y))
  end
end

test_poly_vals_eq(2x + 2(x*y), Polynomial([x, Monomial(1,1)], [2.0, 2.0]))

test_poly_vals_eq((y+x)^2, y^2 + 2(x*y) + x^2)

assert(canonical_form(x + 2x*y + y^2*x + 3.0y^2*x) | string == "1.0 x^1 y^0 + 2.0 x^1 y^1 + 4.0 x^1 y^2")

# polynomial equality
assert((Monomial(0,0)+x)^2 == Monomial(0,0) + 2x + x^2)
assert((y+x)^2 == y^2 + 2(x*y) + x^2)


assert(partial(dim(1), x^2) == 2*x)
assert(partial(dim(1), x^2*y) == 2*x*y)
assert(partial(dim(2), x^2*y) == as_poly(x^2))
assert(partial(dim(2), x^2) == zero_poly(dim(2)))

vm1 = VectorMonomial(x^2*y, dim(2))
div_vm1 = divergence(vm1)
assert(div_vm1 == as_poly(x^2))

vm2 = VectorMonomial(x^2*y^3, dim(1))
div_vm2 = divergence(vm2)
assert(div_vm2 == 2.0 * x*y^3)

assert(Poly.count_monomials_of_degree_eq(pwr(3), dim(2)) == 4) # x^3, x^2 y, x y^2, y^3
assert(count_monomials_of_degree_le(pwr(2), dim(2)) == 6)

mons = monomials_of_degree_le(pwr(2), dim(2))
assert(length(mons) == 6)
assert(mons[1] == Monomial(0,0))
assert(mons[6] == Monomial(2,0))

vmons = vector_monomials_of_degree_le(pwr(2), dim(2))
assert(length(vmons) == 12)
assert(vmons[1] == VectorMonomial(Monomial(0,0), dim(1)))
assert(vmons[2] == VectorMonomial(Monomial(0,0), dim(2)))
assert(vmons[11] == VectorMonomial(Monomial(2,0), dim(1)))
assert(vmons[12] == VectorMonomial(Monomial(2,0), dim(2)))

# dot product of vector monomials
assert(dot(vm1,vm2).coefs == [zeroR])

vm3 = VectorMonomial(x * y, dim(2))
assert(dot(vm1,vm3) == as_poly(x^3 * y^2))

# component extraction of vector monomials
assert(vm1[1] == zero_poly(dim(2)))
assert(vm1[2] == 1*x^2 * y)

# integration of monomials
assert(integrate_on_llo_rect(x*y, 1.0, 1.0) == .25)
assert(integrate_on_llo_rect(x*y, 2.0, 3.0) == 9.)

# integration of polynomials
assert(integrate_on_llo_rect(2(x*y), 1.0, 1.0) == 0.5)
assert(integrate_on_llo_rect(2(x*y) + 4(x*y), 1.0, 1.0) == 1.5)

# 3D polynomial integration

x1 = Monomial(1,0,0)
x2 = Monomial(0,1,0)
x3 = Monomial(0,0,1)

assert(integrate_on_llo_rect(x1*x2*x3, 1.0,2.0,3.0) == 4.5)
assert(integrate_on_llo_rect(x1*x2^2*x3, 1.0,2.0,3.0) == 6.)

# dimension reduction by fixing an input value in one dimension
assert(reduce_dim_by_fixing(dim(2), 2.0, x^2 + -2x*y + y^2) == let x = Monomial(1); x^2 + -4x + 4 end)
assert(reduce_dim_by_fixing(dim(1), 2.0, x^2 + -2x*y + y^2) == let x = Monomial(1); -4x + 4.0 + x^2 end)

# unary minus
assert(-x == -1.0 * x)
assert(-(x^2 + -2x) == 2x + -1x^2)

# polynomial vector operations
pv1 = 1 * vm1
pv2 = 1 * vm2
assert(pv1[2] == as_poly(x^2 * y))
assert(pv1 + pv2 == PolynomialVector([1.0 * x^2 * y^3, 1.0 * x^2 * y]))
assert(2*pv1 + 3*pv2 == PolynomialVector([3.0 * x^2 * y^3, 2.0 * x^2 * y]))

# dot of vector polynomials
assert(dot(pv1,pv1) == 1.0 * x^4 * y^2)

# Test closeness functions
assert(coefs_closer_than(10e-3, 2.0000001x*y + -0.99999999x^2 + 0.00003y^2, 2.0000008x*y + -1.000000001x^2))
assert(!coefs_closer_than(10e-3, 2.0000001x*y + -0.99999999x^2, 2.1000008x*y + -1.000000001x^2))
assert(drop_coefs_lt(10e-3, 0.00000003x*y + 1.2x) == 1.2x)
