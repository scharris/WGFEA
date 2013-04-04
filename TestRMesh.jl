using Test
using RMesh

using Common
import Mesh, Mesh.NBSideInclusions, Mesh.fe_num, Mesh.nb_side_num
import Poly.Monomial, Poly.VectorMonomial

mesh_min_coords = [1.0, 2.0, 3.0]
mesh_max_coords = [4.0, 6.0, 8.0]
mesh_ldims = [mesh_coord(3), mesh_coord(4), mesh_coord(5)]
rmesh3x4x5 = RectMesh(mesh_min_coords, mesh_max_coords, mesh_ldims)

@test rmesh3x4x5.space_dim == 3
@test rmesh3x4x5.min_bounds == mesh_min_coords
@test rmesh3x4x5.max_bounds == mesh_max_coords
@test rmesh3x4x5.mesh_ldims == mesh_ldims
@test rmesh3x4x5.fe_dims == [1.,1.,1.]
@test Mesh.one_mon(rmesh3x4x5) == Monomial(0,0,0)

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
@test Mesh.num_side_faces_per_fe(rmesh3x4x5) == rmesh3x4x5.num_side_faces_per_fe == 6


# Test mesh coordinates and boundary side determination.

@test RMesh.fe_mesh_coord(dim(1), fe_num(1), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(1), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(1), rmesh3x4x5) == 1

@test  Mesh.is_boundary_side(fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(2), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(2), fe_num(2), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(2), rmesh3x4x5) == 1

@test RMesh.fe_mesh_coord(dim(1), fe_num(3), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(3), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(3), rmesh3x4x5) == 1

@test !Mesh.is_boundary_side(fe_num(3), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(3), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(3), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(3), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(3), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(3), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

@test RMesh.fe_mesh_coord(dim(1), fe_num(4), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(4), rmesh3x4x5) == 2
@test RMesh.fe_mesh_coord(dim(3), fe_num(4), rmesh3x4x5) == 1

@test  Mesh.is_boundary_side(fe_num(4), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(4), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(4), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

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

@test !Mesh.is_boundary_side(fe_num(12), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(12), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(12), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(12), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(12), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(12), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)


# first element of second stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(13), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(2), fe_num(13), rmesh3x4x5) == 1
@test RMesh.fe_mesh_coord(dim(3), fe_num(13), rmesh3x4x5) == 2

@test  Mesh.is_boundary_side(fe_num(13), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(13), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(13), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(13), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(13), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(13), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

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

@test  Mesh.is_boundary_side(fe_num(1+4*12), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1+4*12), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(1+4*12), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1+4*12), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(1+4*12), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(1+4*12), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

# last element of fifth and last stack
@test RMesh.fe_mesh_coord(dim(1), fe_num(5*12), rmesh3x4x5) == 3
@test RMesh.fe_mesh_coord(dim(2), fe_num(5*12), rmesh3x4x5) == 4
@test RMesh.fe_mesh_coord(dim(3), fe_num(5*12), rmesh3x4x5) == 5

@test !Mesh.is_boundary_side(fe_num(5*12), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(5*12), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(5*12), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(5*12), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5)
@test !Mesh.is_boundary_side(fe_num(5*12), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5)
@test  Mesh.is_boundary_side(fe_num(5*12), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5)

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
      NBSideInclusions(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(2), RMesh.lesser_side_face_perp_to_axis(dim(1)))

sgeom = RMesh.nb_side_geom(nb_side_num(2), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(1), mesh_coord(1)]


@test Mesh.fe_inclusions_of_nb_side(nb_side_num(2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(2), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(3), RMesh.lesser_side_face_perp_to_axis(dim(1)))

# second row (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(3), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3), rmesh3x4x5) ==
      NBSideInclusions(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(5), RMesh.lesser_side_face_perp_to_axis(dim(1)))

# last side of first stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(8), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(4), mesh_coord(1)]

# last axis-1 perpendicular side in first stack
@test Mesh.fe_inclusions_of_nb_side(nb_side_num(8), rmesh3x4x5) ==
      NBSideInclusions(fe_num(11), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(12), RMesh.lesser_side_face_perp_to_axis(dim(1)))

# first side of second stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(9), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(9), rmesh3x4x5) ==
      NBSideInclusions(fe_num(13), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(14), RMesh.lesser_side_face_perp_to_axis(dim(1)))

# first side of last stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(1+2*4*4), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(5)]

# last side of last stack (axis 1)
sgeom = RMesh.nb_side_geom(nb_side_num(2*4*5), rmesh3x4x5)
@test sgeom.perp_axis == 1
@test sgeom.mesh_coords == [mesh_coord(2), mesh_coord(4), mesh_coord(5)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(2*4*5), rmesh3x4x5) ==
      NBSideInclusions(fe_num(3*4*5-1), RMesh.greater_side_face_perp_to_axis(dim(1)),
                       fe_num(3*4*5), RMesh.lesser_side_face_perp_to_axis(dim(1)))


# Test the non-boundary sides perpendicular to axis 2, which form a 3 x 3 x 5 mesh.

last_axis1 = 2*4*5

# first row (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(4), RMesh.lesser_side_face_perp_to_axis(dim(2)))

sgeom = RMesh.nb_side_geom(nb_side_num(3 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(3), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(6), RMesh.lesser_side_face_perp_to_axis(dim(2)))

# second row (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(4 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(4 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(7), RMesh.lesser_side_face_perp_to_axis(dim(2)))


# last side of first stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(9 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(3), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(9 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(9), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(12), RMesh.lesser_side_face_perp_to_axis(dim(2)))

# first side of second stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(10 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(10 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(13), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(16), RMesh.lesser_side_face_perp_to_axis(dim(2)))

# first side of last stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(1+3*3*4 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(5)]

# last side of last stack (axis 2)
sgeom = RMesh.nb_side_geom(nb_side_num(3*3*5 + last_axis1), rmesh3x4x5)
@test sgeom.perp_axis == 2
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(3), mesh_coord(5)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3*3*5 + last_axis1), rmesh3x4x5) ==
      NBSideInclusions(fe_num(3*4*5-3), RMesh.greater_side_face_perp_to_axis(dim(2)),
                       fe_num(3*4*5), RMesh.lesser_side_face_perp_to_axis(dim(2)))


# Test the non-boundary sides perpendicular to axis 3, which form a 3 x 4 x 4 mesh.

last_axis2 = last_axis1 + 3*3*5

# first row (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(1+12), RMesh.lesser_side_face_perp_to_axis(dim(3)))

sgeom = RMesh.nb_side_geom(nb_side_num(3 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(1), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(3), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(3+12), RMesh.lesser_side_face_perp_to_axis(dim(3)))

# second row (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(2), mesh_coord(1)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(4), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(4+12), RMesh.lesser_side_face_perp_to_axis(dim(3)))

# first side of second stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + 3*4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(2)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + 3*4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(1 + 3*4), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(1 + 3*4 + 12), RMesh.lesser_side_face_perp_to_axis(dim(3)))

# first side of last stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(1 + 3*(3*4) + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(1), mesh_coord(1), mesh_coord(4)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(1 + 3*(3*4) + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(1 + 3*(3*4)), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(1 + 3*(3*4) + 12), RMesh.lesser_side_face_perp_to_axis(dim(3)))


# last side of last stack (axis 3)
sgeom = RMesh.nb_side_geom(nb_side_num(3*4*4 + last_axis2), rmesh3x4x5)
@test sgeom.perp_axis == 3
@test sgeom.mesh_coords == [mesh_coord(3), mesh_coord(4), mesh_coord(4)]

@test Mesh.fe_inclusions_of_nb_side(nb_side_num(3*4*4 + last_axis2), rmesh3x4x5) ==
      NBSideInclusions(fe_num(3*4*5-12), RMesh.greater_side_face_perp_to_axis(dim(3)),
                       fe_num(3*4*5), RMesh.lesser_side_face_perp_to_axis(dim(3)))

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

# Integrate constants on faces of reference element.
@test Mesh.integral_on_ref_fe_face(one_mon, Mesh.interior_face, rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(one_mon, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1.
@test Mesh.integral_on_ref_fe_face(3., Mesh.interior_face, rmesh3x4x5) == 3.
@test Mesh.integral_on_ref_fe_face(3., RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 3.

# Integrate monomials on faces of reference element.

@test Mesh.integral_on_ref_fe_face(y^1*z^2, Mesh.interior_face, rmesh3x4x5) == 1/6
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1/6
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1/6
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 1/3
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(y^1*z^2, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1/2

@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, Mesh.interior_face, rmesh3x4x5) == 1/40
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1/10
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 1/8
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1/20

# Integrate polynomial on faces of reference element.
@test Mesh.integral_on_ref_fe_face(2*y^1*z^2 + 3*x^3*y^4*z^1, Mesh.interior_face, rmesh3x4x5) == 2*(1/6) + 3*(1/40)
@test Mesh.integral_on_ref_fe_face(1.2*y^1*z^2 + 3.4*x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1.2*(1/6)
@test Mesh.integral_on_ref_fe_face(1.5*y^1*z^2 + 0.123*x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1.5*(1/6) + 0.123*(1/10)
@test Mesh.integral_on_ref_fe_face(y^1*z^2 + x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(4.5*y^1*z^2 + x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 4.5*(1/3) + 1/8
@test Mesh.integral_on_ref_fe_face(y^1*z^2 + x^3*y^4*z^1, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_face(y^1*z^2 + x^3*y^4*z^1, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1/2 + 1/20

vm1 = VectorMonomial(x^3*y^4*z, dim(1))
vm2 = VectorMonomial(x^3*y^4*z, dim(2))
vm3 = VectorMonomial(x^3*y^4*z, dim(3))

# Integrate vector monomials vs outward normals on the side faces of the reference finite element.

@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 1/10
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm1, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0

@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 1/8
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm2, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0

@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 0
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(vm3, RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == 1/20

# Integrate on lesser faces using monomials which are constant in the corresponding dimension (so the integrals aren't 0).
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(VectorMonomial(y^4*z,dim(1)), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5) == -1/10
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(VectorMonomial(x^3*z,dim(2)), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5) == -1/8
@test Mesh.integral_on_ref_fe_side_vs_outward_normal(VectorMonomial(x^3*y^4,dim(3)), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5) == -1/20

nearly_eq(a::R, b::R) = abs(b - a) < 10e-8

# Test integration of a product of an arbitrary function and an element-local monomial on finite element faces.

# Integrating the product below on fe 1 should be equivalent to integrating the monomial x^3 y^4 z on the reference element interior.
f1(x::Vector{R}) = (x[1] - mesh_min_coords[1])^2 * (x[2] - mesh_min_coords[2])^3
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), Mesh.interior_face, rmesh3x4x5), 1/40)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5), 1/10)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5), 1/8)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f1, x*y*z, fe_num(1), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5), 1/20)

# Integrating the product below on fe 17 should be equivalent to integrating the monomial x^3 y^4 z on the reference element interior.
fe_coords = RMesh.fe_coords(fe_num(17), rmesh3x4x5)
f2(x::Vector{R}) = (x[1] - fe_coords[1])^2 * (x[2] - fe_coords[2])^3
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), Mesh.interior_face, rmesh3x4x5), 1/40)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.lesser_side_face_perp_to_axis(dim(1)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.greater_side_face_perp_to_axis(dim(1)), rmesh3x4x5), 1/10)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.lesser_side_face_perp_to_axis(dim(2)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.greater_side_face_perp_to_axis(dim(2)), rmesh3x4x5), 1/8)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.lesser_side_face_perp_to_axis(dim(3)), rmesh3x4x5), 0.)
@test nearly_eq(Mesh.integral_prod_on_fe_face(f2, x*y*z, fe_num(17), RMesh.greater_side_face_perp_to_axis(dim(3)), rmesh3x4x5), 1/20)


