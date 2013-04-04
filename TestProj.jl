using Test
using Proj

using Common
import Mesh, Mesh.fe_num
import RMesh, RMesh.RectMesh, RMesh.mesh_coord
import WGBasis.WeakFunsPolyBasis
import Poly.Monomial

mesh_mins = [1.0, 2.0, 3.0]
rmesh3x4x5 = RectMesh(mesh_mins, [4.0, 6.0, 8.0], [mesh_coord(3), mesh_coord(4), mesh_coord(5)])
basis = WeakFunsPolyBasis(deg(4), deg(3), rmesh3x4x5)

x = Monomial(1,0,0)
y = Monomial(0,1,0)
z = Monomial(0,0,1)

# We define a polynomial function of degree 3 on the mesh which is a translation of the monomial 2x^2 y
# to finite element 1's least-coordinates corner. This corner is also the origin of the element's local
# coordinate system for interpreting its local monomials and polynomials.  This function should project
# to the local monomial 2x^2 y on the finite element.

g = (x::Vector{R}) -> 2(x[1] - mesh_mins[1])^2 * (x[2] - mesh_mins[2]) - 3(x[3] - mesh_mins[3])
@test Poly.coefs_closer_than(10e-8, 2x^2*y - 3z, project_onto_fe_face(g, fe_num(1), Mesh.interior_face, basis))
@test Poly.coefs_closer_than(10e-8, 2x^2*y - 3z, project_onto_fe_face(g, fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(1)), basis))
@test Poly.coefs_closer_than(10e-8, 2x^2*y - 3z, project_onto_fe_face(g, fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(1)), basis))




