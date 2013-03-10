using Test
using WGrad

using Common
import Poly, Poly.Monomial, Poly.PolynomialVector
import Mesh
import RMesh2, RMesh2.RectMesh2
import WGBasis.WeakFunsPolyBasis


const k = 4
rmesh = RectMesh2((0.,0.), (3.,3.), 3, 3)
basis = WeakFunsPolyBasis(deg(k), dim(2), rmesh)
wgrad_solver = WGradSolver(deg(k-1), dim(2), rmesh)

x = Monomial(1,0)
y = Monomial(0,1)

# Compute the weak gradient of a weak function which equals a single monomial on the entire finite element,
# so that the weak gradient should equal the standard gradient.
function full_grad_mon(mon::Monomial)
  wg = wgrad(mon, Mesh.interior_face, wgrad_solver) +
       wgrad(mon, RMesh2.top_face, wgrad_solver) +
       wgrad(mon, RMesh2.right_face, wgrad_solver) +
       wgrad(mon, RMesh2.bottom_face, wgrad_solver) +
       wgrad(mon, RMesh2.left_face, wgrad_solver)
  Poly.canonical_form(wg)
end

pv1 = full_grad_mon(x*y)
@test Poly.coefs_closer_than(10e-5, pv1, PolynomialVector([1y,1x]))

pv2 = full_grad_mon(x^2*y)
@test Poly.coefs_closer_than(10e-5, pv2, PolynomialVector([2x*y,1x^2]))

pv3 = full_grad_mon(x^2*y^2)
@test Poly.coefs_closer_than(10e-5, pv3, PolynomialVector([2x*y^2,2x^2*y]))
