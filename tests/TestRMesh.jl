using Test

require("../Common")
require("../Poly")
require("../Mesh")
require("../RMesh")

using RMesh
using Common
import Mesh, Mesh.NBSideInclusions, Mesh.fe_num, Mesh.nb_side_num
import Poly.Monomial, Poly.VectorMonomial

nearly_eq(a::Real, b::Real) = abs(b - a) < 10e-8

mesh_min_coords = [1.0, 2.0, 3.0]
mesh_max_coords = [2.0, 3.0, 4.0]
mesh_ldims = [mesh_coord(3), mesh_coord(4), mesh_coord(5)]
rmesh3x4x5 = RectMesh(mesh_min_coords, mesh_max_coords, mesh_ldims)

@test rmesh3x4x5.space_dim == 3
@test rmesh3x4x5.min_bounds == mesh_min_coords
@test rmesh3x4x5.max_bounds == mesh_max_coords
@test rmesh3x4x5.mesh_ldims == mesh_ldims
@test rmesh3x4x5.fe_dims == [1/3, 1/4, 1/5]
@test Mesh.one_mon(rmesh3x4x5) == Monomial(0,0,0)
@test nearly_eq(Mesh.fe_diameter_inv(fe_num(1), rmesh3x4x5), 1/sqrt((1/3)^2 + (1/4)^2 + (1/5)^2))
@test nearly_eq(Mesh.fe_diameter_inv(fe_num(17), rmesh3x4x5), 1/sqrt((1/3)^2 + (1/4)^2 + (1/5)^2))

@test rmesh3x4x5.cumprods_mesh_ldims == [3, 3*4, 3*4*5]

@test length(rmesh3x4x5.cumprods_nb_side_mesh_ldims_by_perp_axis) == 3
@test rmesh3x4x5.cumprods_nb_side_mesh_ldims_by_perp_axis[1] == [3-1, (3-1)*4, (3-1)*4*5]
@test rmesh3x4x5.cumprods_nb_side_mesh_ldims_by_perp_axis[2] == [3, 3*(4-1), 3*(4-1)*5]
@test rmesh3x4x5.cumprods_nb_side_mesh_ldims_by_perp_axis[3] == [3, 3*4, 3*4*(5-1)]

@test rmesh3x4x5.first_nb_side_nums_by_perp_axis[1] == 1
@test rmesh3x4x5.first_nb_side_nums_by_perp_axis[2] == (3-1)*4*5 + 1
@test rmesh3x4x5.first_nb_side_nums_by_perp_axis[3] == (3-1)*4*5 + 3*(4-1)*5 + 1

@test rmesh3x4x5.num_nb_sides == Mesh.num_nb_sides(rmesh3x4x5) == rmesh3x4x5.first_nb_side_nums_by_perp_axis[3] - 1 + 3 * 4 * 4
@test rmesh3x4x5.num_fes == Mesh.num_fes(rmesh3x4x5) == 3*4*5
@test Mesh.max_num_side_faces_per_fe(rmesh3x4x5) == Mesh.num_side_faces_for_fe(fe_num(1), rmesh3x4x5) == 6
@test Mesh.num_oriented_element_shapes(rmesh3x4x5) == 1

left_face = RMesh.lesser_side_face_perp_to_axis(dim(1))
right_face = RMesh.greater_side_face_perp_to_axis(dim(1))
bottom_face = RMesh.lesser_side_face_perp_to_axis(dim(2))
top_face = RMesh.greater_side_face_perp_to_axis(dim(2))
back_face = RMesh.lesser_side_face_perp_to_axis(dim(3))
front_face = RMesh.greater_side_face_perp_to_axis(dim(3))


# Test mesh coordinates and boundary side determination.

@test RMesh.fe_mesh_coord(dim(1), fe_num(1), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(1), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(1), rmesh3x4x5) == 1

@test  Mesh.Mesh.is_boundary_side(fe_num(1), left_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1), right_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(1), bottom_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1), top_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(1), back_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1), front_face, rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(2), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(2), fe_num(2), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(2), rmesh3x4x5) == 1

@test RMesh.fe_mesh_coord(dim(1), fe_num(3), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(3), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(3), rmesh3x4x5) == 1

@test !Mesh.Mesh.is_boundary_side(fe_num(3), left_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(3), right_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(3), bottom_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(3), top_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(3), back_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(3), front_face, rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(4), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(4), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(3), fe_num(4), rmesh3x4x5) == 1

@test  Mesh.Mesh.is_boundary_side(fe_num(4), left_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(4), right_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(4), bottom_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(4), top_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(4), back_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(4), front_face, rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(6), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(6), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(3), fe_num(6), rmesh3x4x5) == 1

@test RMesh.fe_mesh_coord(dim(1), fe_num(7), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(7), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(3), fe_num(7), rmesh3x4x5) == 1

@test RMesh.fe_mesh_coord(dim(1), fe_num(10), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(10), rmesh3x4x5) == 4
@test RMesh.fe_mesh_coord(dim(3), fe_num(10), rmesh3x4x5) == 1

@test RMesh.fe_mesh_coord(dim(1), fe_num(12), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(12), rmesh3x4x5) == 4
@test RMesh.fe_mesh_coord(dim(3), fe_num(12), rmesh3x4x5) == 1

@test !Mesh.Mesh.is_boundary_side(fe_num(12), left_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(12), right_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(12), bottom_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(12), top_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(12), back_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(12), front_face, rmesh3x4x5)


# first element of second stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(13), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(13), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(13), rmesh3x4x5) == 2

@test  Mesh.Mesh.is_boundary_side(fe_num(13), left_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(13), right_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(13), bottom_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(13), top_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(13), back_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(13), front_face, rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(15), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(15), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(15), rmesh3x4x5) == 2

@test RMesh.fe_mesh_coord(dim(1), fe_num(16), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(16), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(3), fe_num(16), rmesh3x4x5) == 2

# last element of second stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(13+12-1), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(13+12-1), rmesh3x4x5) == 4
@test RMesh.fe_mesh_coord(dim(3), fe_num(13+12-1), rmesh3x4x5) == 2

# first element of third stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(13+12), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(13+12), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(13+12), rmesh3x4x5) == 3

# first element of fifth and last stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(1+4*12), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(1+4*12), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(1+4*12), rmesh3x4x5) == 5

@test  Mesh.Mesh.is_boundary_side(fe_num(1+4*12), left_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1+4*12), right_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(1+4*12), bottom_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1+4*12), top_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(1+4*12), back_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(1+4*12), front_face, rmesh3x4x5)

# last element of fifth and last stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(5*12), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(5*12), rmesh3x4x5) == 4
@test RMesh.fe_mesh_coord(dim(3), fe_num(5*12), rmesh3x4x5) == 5

@test !Mesh.Mesh.is_boundary_side(fe_num(5*12), left_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(5*12), right_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(5*12), bottom_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(5*12), top_face, rmesh3x4x5)
@test !Mesh.Mesh.is_boundary_side(fe_num(5*12), back_face, rmesh3x4x5)
@test  Mesh.Mesh.is_boundary_side(fe_num(5*12), front_face, rmesh3x4x5)

# out of range
@test_fails RMesh.fe_mesh_coord(dim(1), fe_num(5*12+1), rmesh3x4x5) > 0
@test_fails RMesh.fe_mesh_coord(dim(2), fe_num(5*12+1), rmesh3x4x5) > 0
@test_fails RMesh.fe_mesh_coord(dim(3), fe_num(5*12+1), rmesh3x4x5) > 0


# test non-boundary side coordinates

# Test the non-boundary sides perpendicular to axis 1.
# The mesh for non-boundary sides perpendicular to axis 1 has dimensions 2 x 4 x 5.

# first row (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(1), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(1),
                       fe_num(1), right_face,
                       fe_num(2), left_face)

sgeom = RMesh.nb_side_geom(nb_side_num(2), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(1), mesh_coord(1)]


@test Mesh.fe_inclusions_of_nb_side(nb_side_num(2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(2),
                       fe_num(2), right_face,
                       fe_num(3), left_face)

# second row (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(3), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(3),
                       fe_num(4), right_face,
                       fe_num(5), left_face)

# last side of first stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(8), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(4), mesh_coord(1)]

# last axis-1 perpendicular side in first stack
@test Mesh.fe_inclusions_of_nb_side(nb_side_num(8), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(8),
                       fe_num(11), right_face,
                       fe_num(12), left_face)

# first side of second stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(9), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(9), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(9),
                       fe_num(13), right_face,
                       fe_num(14), left_face)

# first side of last stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(1+2*4*4), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(5)]

# last side of last stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(2*4*5), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(4), mesh_coord(5)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(2*4*5), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(2*4*5),
                       fe_num(3*4*5-1), right_face,
                       fe_num(3*4*5), left_face)


# Test the non-boundary sides perpendicular to axis 2, which form a 3 x 3 x 5 mesh.

last_axis1 = 2*4*5

# first row (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(1 + last_axis1),
                       fe_num(1), top_face,
                       fe_num(4), bottom_face)

sgeom = RMesh.nb_side_geom(nb_side_num(3 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(3 + last_axis1),
                       fe_num(3), top_face,
                       fe_num(6), bottom_face)

# second row (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(4 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(4 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(4 + last_axis1),
                       fe_num(4), top_face,
                       fe_num(7), bottom_face)


# last side of first stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(9 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(3), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(9 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(9 + last_axis1),
                       fe_num(9), top_face,
                       fe_num(12), bottom_face)

# first side of second stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(10 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(10 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(10 + last_axis1),
                       fe_num(13), top_face,
                       fe_num(16), bottom_face)

# first side of last stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(1+3*3*4 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(5)]

# last side of last stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(3*3*5 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(3), mesh_coord(5)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3*3*5 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(3*3*5 + last_axis1),
                       fe_num(3*4*5-3), top_face,
                       fe_num(3*4*5), bottom_face)


# Test the non-boundary sides perpendicular to axis 3, which form a 3 x 4 x 4 mesh.

last_axis2 = last_axis1 + 3*3*5

# first row (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(1 + last_axis2),
                       fe_num(1), front_face,
                       fe_num(1+12), back_face)

sgeom = RMesh.nb_side_geom(nb_side_num(3 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(3 + last_axis2),fe_num(3), front_face,
                       fe_num(3+12), back_face)

# second row (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(4 + last_axis2),
                       fe_num(4), front_face,
                       fe_num(4+12), back_face)

# first side of second stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + 3*4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + 3*4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(1 + 3*4 + last_axis2),
                       fe_num(1 + 3*4), front_face,
                       fe_num(1 + 3*4 + 12), back_face)

# first side of last stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + 3*(3*4) + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(4)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + 3*(3*4) + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(1 + 3*(3*4) + last_axis2),
                       fe_num(1 + 3*(3*4)), front_face,
                       fe_num(1 + 3*(3*4) + 12), back_face)


# last side of last stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(3*4*4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(4), mesh_coord(4)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3*4*4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(nb_side_num(3*4*4 + last_axis2),
                       fe_num(3*4*5-12), front_face,
                       fe_num(3*4*5), back_face)

# side number out of range
@test_fails RMesh.nb_side_geom(nb_side_num(1 + 3*4*4 + last_axis2), rmesh3x4x5)


# Test conversion of fe coords to fe number.
@test RMesh.fe_with_mesh_coords([mesh_coord(1), mesh_coord(1), mesh_coord(1)], rmesh3x4x5) == 1
@test RMesh.fe_with_mesh_coords([mesh_coord(3), mesh_coord(1), mesh_coord(1)], rmesh3x4x5) == 3
@test RMesh.fe_with_mesh_coords([mesh_coord(1), mesh_coord(2), mesh_coord(1)], rmesh3x4x5) == 4
@test RMesh.fe_with_mesh_coords([mesh_coord(2), mesh_coord(2), mesh_coord(1)], rmesh3x4x5) == 5
@test RMesh.fe_with_mesh_coords([mesh_coord(1), mesh_coord(1), mesh_coord(2)], rmesh3x4x5) == 13
@test RMesh.fe_with_mesh_coords([mesh_coord(3), mesh_coord(4), mesh_coord(5)], rmesh3x4x5) == 12*5

# Test integrals on finite element faces.

x = Monomial(1,0,0)
y = Monomial(0,1,0)
z = Monomial(0,0,1)
one_mon = RMesh.one_mon(rmesh3x4x5)
VM = VectorMonomial # alias for brevity

rect_oshape = Mesh.oshape(1)

# Integrate constants on faces of reference element.
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  1/3 * 1/4 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, left_face, rmesh3x4x5),
  1/4 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, right_face, rmesh3x4x5),
  1/4 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, bottom_face, rmesh3x4x5),
  1/3 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, top_face, rmesh3x4x5),
  1/3 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, back_face, rmesh3x4x5),
  1/3 * 1/4
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(one_mon, rect_oshape, front_face, rmesh3x4x5),
  1/3 * 1/4
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(3., rect_oshape, Mesh.interior_face, rmesh3x4x5),
  3 * 1/3 * 1/4 * 1/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(3., rect_oshape, front_face, rmesh3x4x5),
  3 * 1/3 * 1/4
)

# Integrate monomials on faces of reference element.

@test nearly_eq(
  Mesh.integral_face_rel_on_face(x*y*z^2, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  (1/3)^2/2 * (1/4)^2/2 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y*z^2, rect_oshape, left_face, rmesh3x4x5),
  (1/4)^2/2 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y*z^2, rect_oshape, right_face, rmesh3x4x5),
  (1/4)^2/2 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x*z^2, rect_oshape, bottom_face, rmesh3x4x5),
  (1/3)^2/2 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y*z^2, rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y*z^2, rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^2*y, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^3/3 * (1/4)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^2*y, rect_oshape, front_face, rmesh3x4x5),
  (1/3)^3/3 * (1/4)^2/2
)

@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^3*y^4*z^1, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y^4*z^1, rect_oshape, left_face, rmesh3x4x5),
  (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y^4*z^1, rect_oshape, right_face, rmesh3x4x5),
  (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^3*z^1, rect_oshape, bottom_face, rmesh3x4x5),
  (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^3*z^1, rect_oshape, top_face, rmesh3x4x5),
  (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^3*y^4, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^3*y^4, rect_oshape, front_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)

# Integrate polynomial on faces of reference element.
@test nearly_eq(
  Mesh.integral_face_rel_on_face(2*y^1*z^2 + 3*x^3*y^4*z^1, rect_oshape, Mesh.interior_face, rmesh3x4x5),
  2(1/3 * (1/4)^2/2 * (1/5)^3/3) + 3((1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(1.2*y^1*z^2 + 3.4*x^3*y^4*z^1, rect_oshape, left_face, rmesh3x4x5),
  1.2((1/4)^2/2 * (1/5)^3/3)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(1.5*y^1*z^2 + 0.123*y^4*z^1, rect_oshape, right_face, rmesh3x4x5),
  1.5((1/4)^2/2 * (1/5)^3/3) + 0.123((1/4)^5/5 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(y^1*z^2 + x^3*y^4*z^1, rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(4.5*x^1*z^2 + 2.3x^3*z^1, rect_oshape, top_face, rmesh3x4x5),
  4.5((1/3)^2/2 * (1/5)^3/3) + 2.3((1/3)^4/4 * (1/5)^2/2)
)
@test nearly_eq(
  Mesh.integral_face_rel_on_face(x^1*y^2 + x^3*y^4, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^2/2 * (1/4)^3/3 + (1/3)^4/4 * (1/4)^5/5
)


# Integrate vector monomials vs outward normals on the side faces of the reference finite element.

@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  (1/3)^3 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(z, VM(y^4*z, dim(1)), rect_oshape, left_face, rmesh3x4x5),
  -(1/4)^5/5 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(y^2*z,   VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  (1/3)^3 * (1/4)^7/7 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x*y^2*z, VM(x^3*y^4*z, dim(1)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(1)), rect_oshape, front_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  (1/4)^4 * (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x^2*z,   VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  (1/4)^4 * (1/3)^6/6 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x^2*y*z, VM(x^3*y^4*z, dim(2)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(2)), rect_oshape, front_face, rmesh3x4x5),
  0
)

@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^4/4 * (1/4)^5/5
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x^2*y,   VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^6/6 * (1/4)^6/6
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x^2*y*z, VM(x^3*y^4*z, dim(3)), rect_oshape, front_face, rmesh3x4x5),
  0
)

# Integrate on lesser faces using monomials which are constant in the corresponding dimension (so the integrals aren't 0).
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(y^4*z,   dim(1)), rect_oshape, left_face, rmesh3x4x5),
  -(1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*z,   dim(2)), rect_oshape, bottom_face, rmesh3x4x5),
  -(1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4, dim(3)), rect_oshape, back_face, rmesh3x4x5),
  -(1/3)^4/4 * (1/4)^5/5
)


# Test integration of a product of an arbitrary function and an element-local monomial on finite element faces.

f1(x::Vector{R}) = (x[1] - mesh_min_coords[1])^2 * (x[2] - mesh_min_coords[2])^3
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y*z, fe_num(1), Mesh.interior_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, y*z, fe_num(1), left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, y*z, fe_num(1), right_face, rmesh3x4x5),
  (1/3)^2 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, x*z, fe_num(1), bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, x*z, fe_num(1), top_face, rmesh3x4x5),
  (1/4)^3 * (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), back_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f1, x*y, fe_num(1), front_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)

fe17_coords = RMesh.fe_coords(fe_num(17), rmesh3x4x5)
f2(x::Vector{R}) = (x[1] - fe17_coords[1])^2 * (x[2] - fe17_coords[2])^3
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y*z, fe_num(17), Mesh.interior_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, y*z, fe_num(17), left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, y*z, fe_num(17), right_face, rmesh3x4x5),
  (1/3)^2 * (1/4)^5/5 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, x*z, fe_num(17), bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, x*z, fe_num(17), top_face, rmesh3x4x5),
  (1/4)^3 * (1/3)^4/4 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(17), back_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)
@test nearly_eq(
  Mesh.integral_global_x_face_rel_on_fe_face(f2, x*y, fe_num(17), front_face, rmesh3x4x5),
  (1/3)^4/4 * (1/4)^5/5
)


@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y, y, rect_oshape, right_face, rmesh3x4x5),
  (1/3) * (1/4)^3/3 * (1/5)
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y*z, y*z, rect_oshape, right_face, rmesh3x4x5),
  (1/3) * (1/4)^3/3 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y*z, x*y*z, rect_oshape, right_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y, y, rect_oshape, left_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y^2*z, y*z, rect_oshape, left_face, rmesh3x4x5),
  (1/4)^4/4 * (1/5)^3/3
)

@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y*z, y, rect_oshape, top_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y*z, x, rect_oshape, top_face, rmesh3x4x5),
  (1/4) * (1/3)^2/2 * (1/5)^2/2
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y*z, x*z, rect_oshape, top_face, rmesh3x4x5),
  (1/4) * (1/3)^3/3 * (1/5)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y, x*z, rect_oshape, bottom_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*z, z, rect_oshape, bottom_face, rmesh3x4x5),
  (1/3)^2/2 * (1/5)^3/3
)

@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y*x, z, rect_oshape, front_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(z*x, y, rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^2/2 * (1/4)^2/2
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(z*x*y, x*y, rect_oshape, front_face, rmesh3x4x5),
  (1/5) * (1/3)^3/3 * (1/4)^3/3
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(y*z, y*x, rect_oshape, back_face, rmesh3x4x5),
  0
)
@test nearly_eq(
  Mesh.integral_fe_rel_x_side_rel_on_side(x*y, x*y^2, rect_oshape, back_face, rmesh3x4x5),
  (1/3)^3/3 * (1/4)^4/4
)

# 10 rows x 20 cols mesh, with vertexes at each integer pair

rmesh20x10 = RectMesh([0.,0.], [20.,10.], [mesh_coord(20), mesh_coord(10)])

@test RMesh.fe_coords(fe_num(1), rmesh20x10) == [0., 0.]
@test RMesh.fe_coords(fe_num(2), rmesh20x10) == [1., 0.]
@test RMesh.fe_coords(fe_num(20), rmesh20x10) == [19., 0.]
@test RMesh.fe_coords(fe_num(21), rmesh20x10) == [0., 1.]
@test RMesh.fe_coords(fe_num(200), rmesh20x10) == [19., 9.]

@test Mesh.num_fes(rmesh20x10) == 200
@test Mesh.num_nb_sides(rmesh20x10) == 370

# test boundary sides

@test !Mesh.is_boundary_side(fe_num(1), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(1), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(1), bottom_face,  rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(1), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fe_num(8), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(8), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(8), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(8), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fe_num(20), top_face, rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(20), right_face, rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(20), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(20), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fe_num(21), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(21), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(21), bottom_face,  rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(21), left_face,  rmesh20x10)

@test !Mesh.is_boundary_side(fe_num(22), top_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(22), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(22), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(22), left_face,  rmesh20x10)

@test  Mesh.is_boundary_side(fe_num(200), top_face, rmesh20x10)
@test  Mesh.is_boundary_side(fe_num(200), right_face, rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(200), bottom_face,  rmesh20x10)
@test !Mesh.is_boundary_side(fe_num(200), left_face,  rmesh20x10)


rmesh3x2 = RectMesh([0.,0.], [3.,2.], [mesh_coord(3), mesh_coord(2)])

@test RMesh.perp_axis_for_nb_side(nb_side_num(1), rmesh3x2) == dim(1)
@test RMesh.perp_axis_for_nb_side(nb_side_num(4), rmesh3x2) == dim(1)
@test RMesh.perp_axis_for_nb_side(nb_side_num(5), rmesh3x2) == dim(2)
@test RMesh.perp_axis_for_nb_side(nb_side_num(7), rmesh3x2) == dim(2)
@test_fails RMesh.perp_axis_for_nb_side(nb_side_num(8), rmesh3x2)

# Test side inclusions
# fe vertical sides
incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(1), rmesh3x2)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(2), rmesh3x2)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(3)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(3), rmesh3x2)
@test incls.fe1 == fe_num(4)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == left_face

incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(4), rmesh3x2)
@test incls.fe1 == fe_num(5)
@test incls.face_in_fe1 == right_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == left_face

# fe horizontal sides
incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(5), rmesh3x2)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == bottom_face

incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(6), rmesh3x2)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == bottom_face

incls = Mesh.fe_inclusions_of_nb_side(nb_side_num(7), rmesh3x2)
@test incls.fe1 == fe_num(3)
@test incls.face_in_fe1 == top_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == bottom_face


# Test integrals on finite element faces.

x = Monomial(1,0)
y = Monomial(0,1)
one_mon = RMesh.one_mon(rmesh3x2)

# Integrate constants on faces of reference element.
@test Mesh.integral_face_rel_on_face(one_mon, rect_oshape, Mesh.interior_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, rect_oshape, left_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, rect_oshape, right_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, rect_oshape, bottom_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_face(one_mon, rect_oshape, top_face, rmesh3x2) == 1.
@test Mesh.integral_face_rel_on_face(3., rect_oshape, Mesh.interior_face, rmesh3x2) == 3.
@test Mesh.integral_face_rel_on_face(3., rect_oshape, top_face, rmesh3x2) == 3.

# Integrate monomials on faces of reference element.

@test Mesh.integral_face_rel_on_face(x*y^2, rect_oshape, Mesh.interior_face, rmesh3x2) == 1/6
@test Mesh.integral_face_rel_on_face(y^2, rect_oshape, left_face, rmesh3x2) == 1/3
@test Mesh.integral_face_rel_on_face(y^2, rect_oshape, right_face, rmesh3x2) == 1/3
@test Mesh.integral_face_rel_on_face(x*y^2, rect_oshape, right_face, rmesh3x2) == 0
@test Mesh.integral_face_rel_on_face(x, rect_oshape, bottom_face, rmesh3x2) == 1/2
@test Mesh.integral_face_rel_on_face(x, rect_oshape, top_face, rmesh3x2) == 1/2
@test Mesh.integral_face_rel_on_face(x*y, rect_oshape, top_face, rmesh3x2) == 0

@test Mesh.integral_face_rel_on_face(x^3*y^4, rect_oshape, Mesh.interior_face, rmesh3x2) == 1/20
@test Mesh.integral_face_rel_on_face(y^4, rect_oshape, left_face, rmesh3x2) == 1/5
@test Mesh.integral_face_rel_on_face(y^4, rect_oshape, right_face, rmesh3x2) == 1/5
@test Mesh.integral_face_rel_on_face(x*y^4, rect_oshape, right_face, rmesh3x2) == 0
@test Mesh.integral_face_rel_on_face(x^3, rect_oshape, bottom_face, rmesh3x2) == 1/4
@test Mesh.integral_face_rel_on_face(x^3, rect_oshape, top_face, rmesh3x2) == 1/4
@test Mesh.integral_face_rel_on_face(x^3*y, rect_oshape, top_face, rmesh3x2) == 0

## Integrate polynomial on faces of reference element.
@test nearly_eq(Mesh.integral_face_rel_on_face(2*x*y^2 + 3*x^3*y^4, rect_oshape, Mesh.interior_face, rmesh3x2), 2*(1/6) + 3*(1/20))
@test nearly_eq(Mesh.integral_face_rel_on_face(1.2*y^2 + 3.4*y^4, rect_oshape, left_face, rmesh3x2), 1.2/3 + 3.4/5)
@test nearly_eq(Mesh.integral_face_rel_on_face(1.5*y^2 + 0.123*y^4, rect_oshape, right_face, rmesh3x2), 1.5*(1/3) + 0.123*(1/5))
@test nearly_eq(Mesh.integral_face_rel_on_face(1.2x + 2x^3, rect_oshape, bottom_face, rmesh3x2), 1.2/2 + 2/4)
@test nearly_eq(Mesh.integral_face_rel_on_face(4.5*x*y^2 + 23.2x^3, rect_oshape, top_face, rmesh3x2), 23.2/4)

# Integrate vector monomials vs outward normals on the side faces of the reference finite element.
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(one_mon, VM(x^3*y^4, dim(1)), rect_oshape, left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(y, VM(x^3*y^4, dim(1)), rect_oshape, right_face, rmesh3x2), 1/6)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(y, VM(y^2, dim(1)), rect_oshape, left_face, rmesh3x2), -1/4)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x*y, VM(y^2, dim(1)), rect_oshape, right_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(y, VM(x*y^2, dim(1)), rect_oshape, bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(y, VM(y^2, dim(1)), rect_oshape, top_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x, VM(x^2, dim(2)), rect_oshape, bottom_face, rmesh3x2), -1/4)
@test nearly_eq(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(x, VM(x^2, dim(2)), rect_oshape, top_face, rmesh3x2), 1/4)

# Test integration of a product of an arbitrary function and an element-local monomial on finite element faces.

f3(x::Vector{R}) = x[1]^2 * x[2]^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x*y, fe_num(1), Mesh.interior_face, rmesh3x2), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, y, fe_num(1), left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, y, fe_num(1), right_face, rmesh3x2), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x, fe_num(1), bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x, fe_num(1), top_face, rmesh3x2), 1/4)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f3, x*y, fe_num(1), top_face, rmesh3x2), 0)

# Integrating the product below on fe 5 should be equivalent to integrating the monomial x^3 y^4 z on the reference element interior.
fe5_coords = RMesh.fe_coords(fe_num(5), rmesh3x2)
f4(x::Vector{R}) = (x[1] - fe5_coords[1])^2 * (x[2] - fe5_coords[2])^3
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x*y, fe_num(5), Mesh.interior_face, rmesh3x2), 1/20)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, y, fe_num(5), left_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, y, fe_num(5), right_face, rmesh3x2), 1/5)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x, fe_num(5), bottom_face, rmesh3x2), 0)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x, fe_num(5), top_face, rmesh3x2), 1/4)
@test nearly_eq(Mesh.integral_global_x_face_rel_on_fe_face(f4, x*y, fe_num(5), top_face, rmesh3x2), 0)
