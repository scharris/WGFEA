using Test
using Proj

using Common
import Mesh, Mesh.fe_num
import RMesh2
import WGBasis.WeakFunsPolyBasis
import Poly.Monomial

basis = WeakFunsPolyBasis(deg(3), deg(2), RMesh2.RectMesh2((0.,0.), (3.,2.), 2, 3))

x = Monomial(1,0)
y = Monomial(0,1)

# identity projection
g(x) = x^2 * y
pg = project_onto_fe_face(g, fe_num(1), Mesh.interior_face, basis)

