using Test
using WGrad

using Common
import Poly, Poly.Monomial, Poly.PolynomialVector, Poly.Nomial
import Mesh
import RMesh2, RMesh2.RectMesh2
import WGBasis.WeakFunsPolyBasis

const k = 4
rmesh = RectMesh2((0.,0.), (3.,3.), 3, 3)
wgrad_solver = WGradSolver(deg(k-1), rmesh)

x = Monomial(1,0)
y = Monomial(0,1)

# Compute the weak gradient of a weak function which equals a single monomial on the entire finite element,
# so that the weak gradient should equal the standard gradient.
function sum_wgrad_all_faces(v::Nomial)
  wg = wgrad(v, Mesh.interior_face, wgrad_solver) +
       wgrad(v, RMesh2.top_face, wgrad_solver) +
       wgrad(v, RMesh2.right_face, wgrad_solver) +
       wgrad(v, RMesh2.bottom_face, wgrad_solver) +
       wgrad(v, RMesh2.left_face, wgrad_solver)
  Poly.canonical_form(wg)
end


# Test wgrad of monomials.
pv1 = sum_wgrad_all_faces(x*y)
@test Poly.coefs_closer_than(10e-5, pv1, PolynomialVector([1y,1x]))

pv2 = sum_wgrad_all_faces(x^2*y)
@test Poly.coefs_closer_than(10e-5, pv2, PolynomialVector([2x*y,1x^2]))

pv3 = sum_wgrad_all_faces(x^2*y^2)
@test Poly.coefs_closer_than(10e-5, pv3, PolynomialVector([2x*y^2,2x^2*y]))


# Test wgrad of polynomials.

@test Poly.coefs_closer_than(10e-5, sum_wgrad_all_faces(2x*y + x^2*y), PolynomialVector([2y+2x*y, 2x+x^2]))

@test Poly.coefs_closer_than(10e-5, sum_wgrad_all_faces(x*y + 3x^2*y), PolynomialVector([y+6x*y, x+3x^2]))
