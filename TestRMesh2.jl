using Test
import RMesh2, RMesh2.RectMesh2, RMesh2.fe_col, RMesh2.fe_row, RMesh2.fe_coords
using Mesh

rmesh200x400 = RectMesh2((0.,0.), (200.,100.), 200, 400)

# Returns the number of the finite element to the left, or 0 if there is none.
fe_left_of(fe::FENum, mesh::RectMesh2) =
  let col = fe_col(fe, mesh)
    if col == 1 0 else fe - 1 end
  end

# Returns the number of the finite element to the right, or 0 if there is none.
fe_right_of(fe::FENum, mesh::RectMesh2) =
  let col = fe_col(fe, mesh)
    if col == mesh.cols 0 else fe + 1 end
  end

# Returns the number of the finite element above the passed finite element, or 0 if there is none.
fe_above(fe::FENum, mesh::RectMesh2) =
  let row = fe_row(fe, mesh)
    if row == mesh.rows 0 else fe + mesh.cols end
  end

# Returns the number of the finite element below the passed finite element, or 0 if there is none.
fe_below(fe::FENum, mesh::RectMesh2) =
  let row = fe_row(fe, mesh)
    if row == 1 0 else fe - mesh.cols end
  end

@test fe_left_of(fe_num(1), rmesh200x400) == 0
@test fe_left_of(fe_num(2), rmesh200x400) == 1
@test fe_left_of(fe_num(400), rmesh200x400) == 399
@test fe_left_of(fe_num(401), rmesh200x400) == 0
@test fe_left_of(fe_num(199*400+1), rmesh200x400) == 0
@test fe_left_of(fe_num(200*400), rmesh200x400) == 200*400-1

@test fe_right_of(fe_num(1), rmesh200x400) == 2
@test fe_right_of(fe_num(400), rmesh200x400) == 0
@test fe_right_of(fe_num(401), rmesh200x400) == 402
@test fe_right_of(fe_num(199*400+1), rmesh200x400) == 199*400+2
@test fe_right_of(fe_num(200*400), rmesh200x400) == 0

@test fe_above(fe_num(1), rmesh200x400) == 401
@test fe_above(fe_num(401), rmesh200x400) == 801
@test fe_above(fe_num(198*400+1), rmesh200x400) == 199*400+1
@test fe_above(fe_num(200*400), rmesh200x400) == 0

@test fe_below(fe_num(401), rmesh200x400) == 1
@test fe_below(fe_num(400), rmesh200x400) == 0
@test fe_below(fe_num(801), rmesh200x400) == 401
@test fe_below(fe_num(199*400+1), rmesh200x400) == 198*400+1
@test fe_below(fe_num(1), rmesh200x400) == 0


# Test finite element plane coordinates

# 10 rows x 20 cols mesh, with vertexes at each integer pair
rmesh10x20 = RectMesh2((0.,0.), (20.,10.), 10, 20)

@test fe_coords(fe_num(1), rmesh10x20) == [0.0,0.0,1.0,1.0]
@test fe_coords(fe_num(2), rmesh10x20) == [1.0,0.0,2.0,1.0]
@test fe_coords(fe_num(20), rmesh10x20) == [19.0,0.0,20.0,1.0]
@test fe_coords(fe_num(21), rmesh10x20) == [0.0,1.0,1.0,2.0]
@test fe_coords(fe_num(200), rmesh10x20) == [19.0,9.0,20.0,10.0]

@test rmesh10x20.bottom_left_x == 0.
@test rmesh10x20.bottom_left_y == 0.
@test rmesh10x20.top_right_x == 20.
@test rmesh10x20.top_right_y == 10.
@test rmesh10x20.rows == 10
@test rmesh10x20.cols == 20
@test rmesh10x20.fe_width == 1.
@test rmesh10x20.fe_height == 1.
@test rmesh10x20.num_elements == 200
@test rmesh10x20.num_nb_vert_sides == 190
@test rmesh10x20.num_nb_horz_sides == 180
@test rmesh10x20.num_nb_sides == 370
@test rmesh10x20.sidenum_first_nb_horz_side == 191

@test num_fes(rmesh10x20) == 200
@test num_nb_sides(rmesh10x20) == 370
@test fe_for_interior(fe_num(23), rmesh10x20) == 23


rmesh2x3 = RMesh2.RectMesh2((0.,0.), (3.,2.), 2, 3)
@test RMesh2.is_vert_nb_side(side_num(1), rmesh2x3)
@test RMesh2.is_vert_nb_side(side_num(4), rmesh2x3)
@test RMesh2.is_horz_nb_side(side_num(5), rmesh2x3)
@test RMesh2.is_horz_nb_side(side_num(7), rmesh2x3)
@test RMesh2.is_horz_nb_side(side_num(7), rmesh2x3)
@test !RMesh2.is_horz_nb_side(side_num(8), rmesh2x3)

# Test side inclusions
incls = Mesh.NBSideInclusions()

# fe vertical sides
fe_inclusions_of_nb_side!(side_num(1), rmesh2x3, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(side_num(2), rmesh2x3, incls)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(3)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(side_num(3), rmesh2x3, incls)
@test incls.fe1 == fe_num(4)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == RMesh2.left_face

fe_inclusions_of_nb_side!(side_num(4), rmesh2x3, incls)
@test incls.fe1 == fe_num(5)
@test incls.face_in_fe1 == RMesh2.right_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == RMesh2.left_face

# fe horizontal sides
fe_inclusions_of_nb_side!(side_num(5), rmesh2x3, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(4)
@test incls.face_in_fe2 == RMesh2.bottom_face

fe_inclusions_of_nb_side!(side_num(6), rmesh2x3, incls)
@test incls.fe1 == fe_num(2)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(5)
@test incls.face_in_fe2 == RMesh2.bottom_face

fe_inclusions_of_nb_side!(side_num(7), rmesh2x3, incls)
@test incls.fe1 == fe_num(3)
@test incls.face_in_fe1 == RMesh2.top_face
@test incls.fe2 == fe_num(6)
@test incls.face_in_fe2 == RMesh2.bottom_face
