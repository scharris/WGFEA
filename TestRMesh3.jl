using Test
using RMesh3
using Mesh

rmesh200x400x600 = RectMesh3((1.,2.,3.), (201.,402.,603.), 200, 400, 600)

@test rmesh200x400x600.min_x == 1.
@test rmesh200x400x600.min_y == 2.
@test rmesh200x400x600.min_z == 3.

@test rmesh200x400x600.max_x == 201.
@test rmesh200x400x600.max_y == 402.
@test rmesh200x400x600.max_z == 603.

@test rmesh200x400x600.cols == 200
@test rmesh200x400x600.rows == 400
@test rmesh200x400x600.stacks == 600


@test rmesh200x400x600.fe_width  == 1.
@test rmesh200x400x600.fe_height == 1.
@test rmesh200x400x600.fe_depth  == 1.

@test rmesh200x400x600.num_elements == 200 * 400 * 600

@test rmesh200x400x600.num_x_nb_sides == 199 * 400 * 600
@test rmesh200x400x600.num_y_nb_sides == 200 * 399 * 600
@test rmesh200x400x600.num_z_nb_sides == 200 * 400 * 599

@test rmesh200x400x600.num_nb_sides == 199*400*600 + 200*399*600 + 200*400*599

@test rmesh200x400x600.sidenum_first_x_nb_side == 1
@test rmesh200x400x600.sidenum_first_y_nb_side == 199*400*600 + 1
@test rmesh200x400x600.sidenum_first_z_nb_side == 199*400*600 + 200*399*600 + 1

@test rmesh200x400x600.num_fes_per_stack == 200 * 400
@test rmesh200x400x600.num_x_nb_sides_per_stack== 199*400
@test rmesh200x400x600.num_y_nb_sides_per_stack== 200*399
@test rmesh200x400x600.num_z_nb_sides_per_stack== 200*400

@test num_fes(rmesh200x400x600) == 200 * 400 * 600
@test num_nb_sides(rmesh200x400x600) == 199*400*600 + 200*399*600 + 200*400*599

incls = Mesh.NBSideInclusions()

fe_inclusions_of_nb_side!(side_num(1), rmesh200x400x600, incls)
@test incls.fe1 == fe_num(1)
@test incls.face_in_fe1 == RMesh3.x_max_face
@test incls.fe2 == fe_num(2)
@test incls.face_in_fe2 == RMesh3.x_min_face

fe_inclusions_of_nb_side!(side_num(199), rmesh200x400x600, incls)
@test incls.fe1 == fe_num(199)
@test incls.face_in_fe1 == RMesh3.x_max_face
@test incls.fe2 == fe_num(200)
@test incls.face_in_fe2 == RMesh3.x_min_face

fe_inclusions_of_nb_side!(side_num(200), rmesh200x400x600, incls)
@test incls.fe1 == fe_num(201)
@test incls.face_in_fe1 == RMesh3.x_max_face
@test incls.fe2 == fe_num(202)
@test incls.face_in_fe2 == RMesh3.x_min_face


m0 = [1.,2.,3.] # mesh origin

# first x-perpendicular nb side
@test side_info(side_num(1), rmesh200x400x600) ==
      SideInfo(side_num(1), 'x',
        FEInfo(fe_num(1), 1,1,1, m0),
        FEInfo(fe_num(2), 2,1,1, [1.,0.,0.]+m0))

# second x-perpendicular nb side
@test side_info(side_num(2), rmesh200x400x600) ==
      SideInfo(side_num(2), 'x',
        FEInfo(fe_num(2), 2,1,1, [1.,0.,0.]+m0),
        FEInfo(fe_num(3), 3,1,1, [2.,0.,0.]+m0))

# last x-perpendicular nb side in first row
@test side_info(side_num(199), rmesh200x400x600) ==
      SideInfo(side_num(199), 'x',
        FEInfo(fe_num(199), 199,1,1, [198.,0.,0.]+m0),
        FEInfo(fe_num(200), 200,1,1, [199.,0.,0.]+m0))

# first x-perpendicular nb side in second row
@test side_info(side_num(200), rmesh200x400x600) ==
      SideInfo(side_num(200), 'x',
        FEInfo(fe_num(201), 1,2,1, [0.,1.,0.]+m0),
        FEInfo(fe_num(202), 2,2,1, [1.,1.,0.]+m0))

# last x-perpendicular nb side in second row
@test side_info(side_num(398), rmesh200x400x600) ==
      SideInfo(side_num(398), 'x',
        FEInfo(fe_num(399), 199,2,1, [198.,1.,0.]+m0),
        FEInfo(fe_num(400), 200,2,1, [199.,1.,0.]+m0))

# last x-perpendicular nb side in the first stack
@test side_info(side_num(199*400), rmesh200x400x600) ==
      SideInfo(side_num(199*400), 'x',
        FEInfo(fe_num(200*400-1), 199,400,1, [198.,399.,0.]+m0),
        FEInfo(fe_num(200*400),   200,400,1, [199.,399.,0.]+m0))

# first x-perpendicular nb side in the second stack
@test side_info(side_num(199*400+1), rmesh200x400x600) ==
      SideInfo(side_num(199*400+1), 'x',
        FEInfo(fe_num(200*400+1), 1,1,2, [0.,0.,1.]+m0),
        FEInfo(fe_num(200*400+2), 2,1,2, [1.,0.,1.]+m0))

# second x-perpendicular nb side in the second stack
@test side_info(side_num(199*400+2), rmesh200x400x600) ==
      SideInfo(side_num(199*400+2), 'x',
        FEInfo(fe_num(200*400+2), 2,1,2, [1.,0.,1.]+m0),
        FEInfo(fe_num(200*400+3), 3,1,2, [2.,0.,1.]+m0))

# last x-perpendicular nb side in the last stack
@test side_info(side_num(199*400*600), rmesh200x400x600) ==
      SideInfo(side_num(199*400*600), 'x',
        FEInfo(fe_num(200*400*600-1), 199,400,600, [198.,399.,599.]+m0),
        FEInfo(fe_num(200*400*600),   200,400,600, [199.,399.,599.]+m0))


yside(i::Integer) = side_num(199*400*600) + i

# first y-perpendicular nb side in the first stack
@test side_info(yside(1), rmesh200x400x600) ==
      SideInfo(yside(1), 'y',
        FEInfo(fe_num(1),   1,1,1, [0.,0.,0.]+m0),
        FEInfo(fe_num(201), 1,2,1, [0.,1.,0.]+m0))

# second y-perpendicular nb side in the first stack
@test side_info(yside(2), rmesh200x400x600) ==
      SideInfo(yside(2), 'y',
        FEInfo(fe_num(2),   2,1,1, [1.,0.,0.]+m0),
        FEInfo(fe_num(202), 2,2,1, [1.,1.,0.]+m0))

# last y-perpendicular nb side in the first row of the first stack
@test side_info(yside(200), rmesh200x400x600) ==
      SideInfo(yside(200), 'y',
        FEInfo(fe_num(200), 200,1,1, [199.,0.,0.]+m0),
        FEInfo(fe_num(400), 200,2,1, [199.,1.,0.]+m0))

# first y-perpendicular nb side in the second row of the first stack
@test side_info(yside(201), rmesh200x400x600) ==
      SideInfo(yside(201), 'y',
        FEInfo(fe_num(201), 1,2,1, [0.,1.,0.]+m0),
        FEInfo(fe_num(401), 1,3,1, [0.,2.,0.]+m0))

# last y-perpendicular nb side in the first stack
@test side_info(yside(200*399), rmesh200x400x600) ==
      SideInfo(yside(200*399), 'y',
        FEInfo(fe_num(200*399), 200,399,1, [199.,398.,0.]+m0),
        FEInfo(fe_num(200*400), 200,400,1, [199.,399.,0.]+m0))

# first y-perpendicular nb side in the second stack
@test side_info(yside(200*399+1), rmesh200x400x600) ==
      SideInfo(yside(200*399+1), 'y',
        FEInfo(fe_num(200*400+1),   1,1,2, [0.,0.,1.]+m0),
        FEInfo(fe_num(200*400+201), 1,2,2, [0.,1.,1.]+m0))

# second y-perpendicular nb side in the second stack
@test side_info(yside(200*399+1), rmesh200x400x600) ==
      SideInfo(yside(200*399+1), 'y',
        FEInfo(fe_num(200*400+1),   1,1,2, [0.,0.,1.]+m0),
        FEInfo(fe_num(200*400+201), 1,2,2, [0.,1.,1.]+m0))

# last y-perpendicular nb side in the second stack
@test side_info(yside(200*399 + 200*399), rmesh200x400x600) ==
      SideInfo(yside(200*399 + 200*399), 'y',
        FEInfo(fe_num(200*400 + 200*399), 200,399,2, [199.,398.,1.]+m0),
        FEInfo(fe_num(200*400 + 200*400), 200,400,2, [199.,399.,1.]+m0))

# last y-perpendicular nb side in the last stack
@test side_info(yside(200*399*600), rmesh200x400x600) ==
      SideInfo(yside(200*399*600), 'y',
        FEInfo(fe_num(200*400*600-200), 200,399,600, [199.,398.,599.]+m0),
        FEInfo(fe_num(200*400*600),     200,400,600, [199.,399.,599.]+m0))

zside(i::Integer) = side_num(199*400*600) + side_num(200*399*600) + i

# first z-perpendicular nb side in the first stack
@test side_info(zside(1), rmesh200x400x600) ==
      SideInfo(zside(1), 'z',
        FEInfo(fe_num(1),         1,1,1, [0.,0.,0.]+m0),
        FEInfo(fe_num(1+200*400), 1,1,2, [0.,0.,1.]+m0))

# second z-perpendicular nb side in the first stack
@test side_info(zside(2), rmesh200x400x600) ==
      SideInfo(zside(2), 'z',
        FEInfo(fe_num(2),         2,1,1, [1.,0.,0.]+m0),
        FEInfo(fe_num(2+200*400), 2,1,2, [1.,0.,1.]+m0))

# last z-perpendicular nb side in the first row of the first stack
@test side_info(zside(200), rmesh200x400x600) ==
      SideInfo(zside(200), 'z',
        FEInfo(fe_num(200),         200,1,1, [199.,0.,0.]+m0),
        FEInfo(fe_num(200+200*400), 200,1,2, [199.,0.,1.]+m0))

# first z-perpendicular nb side in the second row of the first stack
@test side_info(zside(201), rmesh200x400x600) ==
      SideInfo(zside(201), 'z',
        FEInfo(fe_num(201),         1,2,1, [0.,1.,0.]+m0),
        FEInfo(fe_num(201+200*400), 1,2,2, [0.,1.,1.]+m0))

