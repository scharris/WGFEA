
using Poly
using RMesh
using WGBasis
using WGrad

k = 4
rmesh = RectMesh((0.,0.), (3.,3.), mesh_rc(3), mesh_rc(3))
basis = WeakFunsPolyBasis(pwr(k), dim(2), rmesh)
wgrad_solver = WGradSolver(pwr(k-1), dim(2), rmesh)

x = Monomial(1,0)
y = Monomial(0,1)

# Compute the weak gradient of a weak function which equals a single monomial on the entire finite element,
# so that the weak gradient should equal the standard gradient.
function full_grad_mon(mon::Monomial)
  wg = wgrad(mon, RMesh.interior_face, wgrad_solver) +
       wgrad(mon, RMesh.top_face, wgrad_solver) +
       wgrad(mon, RMesh.right_face, wgrad_solver) +
       wgrad(mon, RMesh.bottom_face, wgrad_solver) +
       wgrad(mon, RMesh.left_face, wgrad_solver)
  canonical_form(wg)
end

pv1 = full_grad_mon(x*y)
assert(coefs_closer_than(10e-5, pv1, PolynomialVector([1y,1x])))

pv2 = full_grad_mon(x^2*y)
assert(coefs_closer_than(10e-5, pv2, PolynomialVector([2x*y,1x^2])))

pv3 = full_grad_mon(x^2*y^2)
assert(coefs_closer_than(10e-5, pv3, PolynomialVector([2x*y^2,2x^2*y])))
