using RMesh

rmesh = RectMesh((0.,0.), (200.,100.), fe_coord(400), fe_coord(200))

assert(rmesh.fe_left_of(fe_num(1)) == 0)
assert(rmesh.fe_left_of(fe_num(2)) == 1)
assert(rmesh.fe_left_of(fe_num(400)) == 399)
assert(rmesh.fe_left_of(fe_num(401)) == 0)
assert(rmesh.fe_left_of(fe_num(199*400+1)) == 0)
assert(rmesh.fe_left_of(fe_num(200*400)) == 200*400-1)

assert(rmesh.fe_right_of(fe_num(1)) == 2)
assert(rmesh.fe_right_of(fe_num(400)) == 0)
assert(rmesh.fe_right_of(fe_num(401)) == 402)
assert(rmesh.fe_right_of(fe_num(199*400+1)) == 199*400+2)
assert(rmesh.fe_right_of(fe_num(200*400)) == 0)

assert(rmesh.fe_above(fe_num(1)) == 401)
assert(rmesh.fe_above(fe_num(401)) == 801)
assert(rmesh.fe_above(fe_num(198*400+1)) == 199*400+1)
assert(rmesh.fe_above(fe_num(200*400)) == 0)

assert(rmesh.fe_below(fe_num(401)) == 1)
assert(rmesh.fe_below(fe_num(400)) == 0)
assert(rmesh.fe_below(fe_num(801)) == 401)
assert(rmesh.fe_below(fe_num(199*400+1)) == 198*400+1)
assert(rmesh.fe_below(fe_num(1)) == 0)
