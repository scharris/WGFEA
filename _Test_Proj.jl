using Test
using Proj
using Common
import Poly.Monomial, Poly.Polynomial
import Mesh, Mesh.fe_num
import RMesh, RMesh.RectMesh, RMesh.mesh_coord
import WGBasis.WeakFunsPolyBasis


mesh_mins = [1.0, 2.0, 3.0]
rmesh3x4x5 = RectMesh(mesh_mins, [4.0, 6.0, 8.0], [mesh_coord(3), mesh_coord(4), mesh_coord(5)])
basis = WeakFunsPolyBasis(deg(5), deg(4), rmesh3x4x5)

x = Monomial(1,0,0)
y = Monomial(0,1,0)
z = Monomial(0,0,1)


left_face = RMesh.lesser_side_face_perp_to_axis(dim(1))
right_face = RMesh.greater_side_face_perp_to_axis(dim(1))
bottom_face = RMesh.lesser_side_face_perp_to_axis(dim(2))
top_face = RMesh.greater_side_face_perp_to_axis(dim(2))
back_face = RMesh.lesser_side_face_perp_to_axis(dim(3))
front_face = RMesh.greater_side_face_perp_to_axis(dim(3))

rect_oshape = Mesh.oshape(1)

# Project a global function onto the interior which should match the local polynomial
# (x,y,z) -> 2x^2*y - 3z there, to which it should project.

f = (x::Vector{R}) -> 2(x[1] - mesh_mins[1])^2 * (x[2] - mesh_mins[2]) - 3(x[3] - mesh_mins[3])
@test Poly.coefs_closer_than(10e-7,
  2x^2 * y - 3z,
  Polynomial(basis.interior_mons, project_onto_fe_face(f, fe_num(1), Mesh.interior_face, basis))
)


# constant projection
@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_fe_side(fe_num(1), left_face, basis), project_onto_fe_face(2.3, fe_num(1), left_face, basis)),
  2.3*Mesh.one_mon(rmesh3x4x5)
)

# Project a global function onto the left side which should match the local polynomial
# (x,y,z) -> 2y^3 - 3.2z^2 there, to which it should project.
f = x::Vector{R} -> 2(x[2] - mesh_mins[2])^3 - 3.2(x[3] - mesh_mins[3])^2 +
                    (x[1] - mesh_mins[1])^3 # (this term should vanish on the side)
@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_fe_side(fe_num(1), left_face, basis), project_onto_fe_face(f, fe_num(1), left_face, basis)),
  2y^3 - 3.2z^2
)

# Project a global function onto the right side which should match the local polynomial
# (x,y,z) -> y^2 - 3z there, to which it should project.
f = x::Vector{R} -> (x[2] - mesh_mins[2])^2 - 3(x[3] - mesh_mins[3]) +
                    (x[1] - (mesh_mins[1]+1))^3 # (this term should vanish on the side)
@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_fe_side(fe_num(1), right_face, basis), project_onto_fe_face(f, fe_num(1), right_face, basis)),
  y^2 - 3z
)

# Project a global function onto the bottom side which should match the local polynomial
# (x,y,z) -> x^3 - 3z^3 there, to which it should project.
f = x::Vector{R} -> (x[1] - mesh_mins[1])^3 - 3(x[3] - mesh_mins[3])^3 +
                    (x[2] - mesh_mins[2])^5 # (this term should vanish on the side)
@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_fe_side(fe_num(1), bottom_face, basis), project_onto_fe_face(f, fe_num(1), bottom_face, basis)),
  x^3 - 3z^3
)

# Project a global function onto the top side which should match the local polynomial
# (x,y,z) -> 2x^2 - 3z^3 + 12 there, to which it should project.
f = x::Vector{R} -> 2(x[1] - mesh_mins[1])^2 - 3(x[3] - mesh_mins[3])^3 + 12 +
                    (x[2] - (mesh_mins[2]+1))^5 # (this term should vanish on the side)
@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_fe_side(fe_num(1), top_face, basis), project_onto_fe_face(f, fe_num(1), top_face, basis)),
  2x^2 - 3z^3 + 12
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, right_face, basis),
             project_interior_mon_onto_oshape_side(x * y^2 * z, rect_oshape, right_face, basis)),
  1y^2 * z
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, left_face, basis),
             project_interior_mon_onto_oshape_side(x * y^2 * z, rect_oshape, left_face, basis)),
  0*x
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, left_face, basis),
             project_interior_mon_onto_oshape_side(y^2 * z, rect_oshape, left_face, basis)),
  1y^2 * z
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, top_face, basis),
             project_interior_mon_onto_oshape_side(x * y^2 * z, rect_oshape, top_face, basis)),
  1(x * z)
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, bottom_face, basis),
             project_interior_mon_onto_oshape_side(x * y^2 * z, rect_oshape, bottom_face, basis)),
  0*x
)

@test Poly.coefs_closer_than(10e-7,
  Polynomial(WGBasis.side_mons_for_oshape_side(rect_oshape, bottom_face, basis),
             project_interior_mon_onto_oshape_side(x * z, rect_oshape, bottom_face, basis)),
  1x*z
)

