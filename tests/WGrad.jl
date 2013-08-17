using Base.Test
using WGrad
using Common
import Poly, Poly.Monomial, Poly.PolynomialVector, Poly.Nomial
import Mesh
import RMesh, RMesh.RectMesh, RMesh.mesh_coord

k = 4
rmesh = RectMesh([0.,0.], [3.,3.], [mesh_coord(3), mesh_coord(3)])
rect_oshape = Mesh.oshapenum(1)
wgrad_solver = WGradSolver(deg(k-1), rmesh)

x = Monomial(1,0)
y = Monomial(0,1)
one = Monomial(0,0)
zero = 0*one

left_face = RMesh.lesser_side_face_perp_to_axis(dim(1))
right_face = RMesh.greater_side_face_perp_to_axis(dim(1))
bottom_face = RMesh.lesser_side_face_perp_to_axis(dim(2))
top_face = RMesh.greater_side_face_perp_to_axis(dim(2))

# Compute the weak gradient of a weak function which equals a single monomial on each face of a finite element,
# so that the weak gradient should equal the standard gradient.
function sum_wgrad_all_faces(v::Nomial)
  wg = wgrad(v, rect_oshape, Mesh.interior_face, wgrad_solver) +
       wgrad(v, rect_oshape, top_face, wgrad_solver) +
       wgrad(v, rect_oshape, right_face, wgrad_solver) +
       wgrad(v, rect_oshape, bottom_face, wgrad_solver) +
       wgrad(v, rect_oshape, left_face, wgrad_solver)
  Poly.canonical_form(wg)
end


# Test wgrad of monomials.

@test Poly.coefs_closer_than(10e-5,
  wgrad(x*y, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x, rect_oshape, top_face, wgrad_solver) +
    wgrad(y, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([1y,1x])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2, rect_oshape, top_face, wgrad_solver) +
    wgrad(y, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y,1x^2])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y + 2.3, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + 2.3, rect_oshape, top_face, wgrad_solver) +
    wgrad(y + 2.3, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + 2.3, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + 2.3, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y,1x^2])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y^2 + -2., rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + -2., rect_oshape, top_face, wgrad_solver) +
    wgrad(y^2 + -2., rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + -2., rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + -2., rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y^2,2x^2*y])
)

@test Poly.coefs_closer_than(10e-5,
  wgrad(x^2*y^2 + -2x*y + 1, rect_oshape, Mesh.interior_face, wgrad_solver) +
    wgrad(x^2 + -2x + 1, rect_oshape, top_face, wgrad_solver) +
    wgrad(y^2 + -2y + 1, rect_oshape, right_face, wgrad_solver) +
    wgrad(zero + 1, rect_oshape, bottom_face, wgrad_solver) +
    wgrad(zero + 1, rect_oshape, left_face, wgrad_solver),
  PolynomialVector([2x*y^2 + -2y,2x^2*y + -2x])
)
