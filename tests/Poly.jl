require("Poly.jl")

using Base.Test
using Common
using Poly

x = Monomial(1,0)
y = Monomial(0,1)

@test x * y == Monomial(1,1)
@test x * y^2 == Monomial(1,2)
@test domain_dim(x) == 2

@test value_at(x * y^2, [1., 2.]) == 4.
@test value_at(2x + 2(x*y), [1.0, 2.0]) == 6.


function test_poly_vals_eq(p1::Polynomial, p2::Polynomial)
  for x in [0.0, 1.0, 2.0, 3.0, 4.0], y in [10.0, 11.0, 12.0, 13.0]
    @test value_at(p1, [x, y]) == value_at(p2, [x, y])
  end
end

test_poly_vals_eq(2x + 2(x*y), Polynomial([x, Monomial(1,1)], [2.0, 2.0]))

test_poly_vals_eq((y+x)^2, y^2 + 2(x*y) + x^2)

@test canonical_form(x + 2x*y + y^2*x + 3.0y^2*x) |> string == "1.0 x^1 y^0 + 2.0 x^1 y^1 + 4.0 x^1 y^2"

# polynomial equality
@test (Monomial(0,0)+x)^2 == Monomial(0,0) + 2x + x^2
@test (y+x)^2 == y^2 + 2(x*y) + x^2


@test partial(dim(1), x^2) == 2*x
@test partial(dim(1), x^2*y) == 2*x*y
@test partial(dim(2), x^2*y) == as_poly(x^2)
@test partial(dim(2), x^2) == zero_poly(dim(2))

@test grad(x^2 * y) == PolynomialVector([2.0*x*y, 1.0*x^2])
@test grad(x^2 * y - 3x*y^3) == PolynomialVector([2.0*x*y - 3y^3, 1.0*x^2 - 9x*y^2])

@test antideriv(dim(1), x) == .5x^2
@test antideriv(dim(2), x) == 1.0x*y
@test antideriv(dim(2), y) == .5y^2
@test antideriv(dim(1), y) == 1.0x*y
@test antideriv(dim(1), 1x^2 + 2x) == x^2 + 1/3 * x^3

vm1 = VectorMonomial(x^2*y, dim(2))
div_vm1 = divergence(vm1)
@test div_vm1 == as_poly(x^2)

vm2 = VectorMonomial(x^2*y^3, dim(1))
div_vm2 = divergence(vm2)
@test div_vm2 == 2.0 * x*y^3

#@test Poly.count_mons_of_deg_eq(deg(3), dim(2)) == 4 # x^3, x^2 y, x y^2, y^3
@test count_mons_of_deg_le(deg(2), dim(2)) == 6
@test length(mons_of_deg_le(deg(5), dim(4))) == count_mons_of_deg_le(deg(5), dim(4))
@test length(mons_of_deg_le(deg(6), dim(5))) == count_mons_of_deg_le(deg(6), dim(5))
# If we hold one dimension constant (with exponent 0), we should get the same count as reducing the domain dimension by 1.
@test length(mons_of_deg_le_with_const_dim(deg(3), dim(1), dim(3))) == count_mons_of_deg_le(deg(3), dim(2))

mons = mons_of_deg_le(deg(2), dim(2))
@test length(mons) == 6
@test mons[1] == Monomial(0,0)
@test mons[6] == Monomial(2,0)

vmons = vector_mons_of_deg_le(deg(2), dim(2))
@test length(vmons) == 12
@test vmons[1] == VectorMonomial(Monomial(0,0), dim(1))
@test vmons[2] == VectorMonomial(Monomial(0,0), dim(2))
@test vmons[11] == VectorMonomial(Monomial(2,0), dim(1))
@test vmons[12] == VectorMonomial(Monomial(2,0), dim(2))

# dot product of vector monomials
@test dot(vm1,vm2).coefs == [zeroR]

vm3 = VectorMonomial(x * y, dim(2))
@test dot(vm1,vm3) == as_poly(x^3 * y^2)

# component extraction of vector monomials
@test vm1[1] == zeroR
@test vm1[2] == x^2 * y

# monomial times vector monomial
@test vm2[1] * vm1 == VectorMonomial(x^4 * y^4, dim(2))
@test vm1[2] * vm1 == VectorMonomial(x^4 * y^2, dim(2))

# integration of monomials
@test integral_on_rect_at_origin(x*y, [1.0, 1.0]) == .25
@test integral_on_rect_at_origin(x*y, [2.0, 3.0]) == 9.

# integration of polynomials
@test integral_on_rect_at_origin(2(x*y), [1.0, 1.0]) == 0.5
@test integral_on_rect_at_origin(2(x*y) + 4(x*y), [1.0, 1.0]) == 1.5

# 3D polynomial integration

x1 = Monomial(1,0,0)
x2 = Monomial(0,1,0)
x3 = Monomial(0,0,1)

@test integral_on_rect_at_origin(x1*x2*x3, [1.0,2.0,3.0]) == 4.5
@test integral_on_rect_at_origin(x1*x2^2*x3, [1.0,2.0,3.0]) == 6.

# dimension reduction by fixing an input value in one dimension
@test reduce_dim_by_fixing(dim(2), 2.0, x^2 + -2x*y + y^2) == let x = Monomial(1); x^2 + -4x + 4 end
@test reduce_dim_by_fixing(dim(1), 2.0, x^2 + -2x*y + y^2) == let x = Monomial(1); -4x + 4.0 + x^2 end

t = Monomial(1)
@test reduce_dim_by_subst(dim(2), t^2 + 1, x*y) == t + t^3
@test reduce_dim_by_subst(dim(1), t^2-2t, x^2*y) == (t^2 - 4t + 4)*t^2*t

poly_path = [t + 1, t^2 - 2t]
@test precompose_with_poly_path(x*y^2, poly_path) == (t+1) * (t^2 - 2t)^2

poly_path = [t + 1, t^2 - 2t, t-3.]
@test precompose_with_poly_path(x1*x2^2*x3^3, poly_path) == (t+1) * (t^2 - 2t)^2 * (t-3)^3

# unary minus
@test -x == -1.0 * x
@test -(x^2 + -2x) == 2x + -1x^2

# polynomial vector operations
pv1 = 1 * vm1
pv2 = 1 * vm2
@test pv1[2] == as_poly(x^2 * y)
@test pv1 + pv2 == PolynomialVector([1.0 * x^2 * y^3, 1.0 * x^2 * y])
@test 2*pv1 + 3*pv2 == PolynomialVector([3.0 * x^2 * y^3, 2.0 * x^2 * y])

# dot of vector polynomials
@test dot(pv1,pv1) == 1.0 * x^4 * y^2

# Test closeness functions.
@test coefs_closer_than(10e-3, 2.0000001x*y + -0.99999999x^2 + 0.00003y^2, 2.0000008x*y + -1.000000001x^2)
@test !coefs_closer_than(10e-3, 2.0000001x*y + -0.99999999x^2, 2.1000008x*y + -1.000000001x^2)
@test drop_coefs_lt(10e-3, 0.00000003x*y + 1.2x) == 1.2x


# Test translation.
@test translate(x^2 * y, [2.,3.]) == (x+2)^2 * (y+3)
@test translate(x * y^2, [-2.,3.]) == (x-2) * (y+3)^2
@test translate(x^3 * y^2, [-2.,3.]) == (x-2)^3 * (y+3)^2
