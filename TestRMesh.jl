using RMesh

import Poly.R

rmesh = RectMesh((0.,0.), (200.,100.), mesh_rc(200), mesh_rc(400))

assert(fe_left_of(fe_num(1), rmesh) == 0)
assert(fe_left_of(fe_num(2), rmesh) == 1)
assert(fe_left_of(fe_num(400), rmesh) == 399)
assert(fe_left_of(fe_num(401), rmesh) == 0)
assert(fe_left_of(fe_num(199*400+1), rmesh) == 0)
assert(fe_left_of(fe_num(200*400), rmesh) == 200*400-1)

assert(fe_right_of(fe_num(1), rmesh) == 2)
assert(fe_right_of(fe_num(400), rmesh) == 0)
assert(fe_right_of(fe_num(401), rmesh) == 402)
assert(fe_right_of(fe_num(199*400+1), rmesh) == 199*400+2)
assert(fe_right_of(fe_num(200*400), rmesh) == 0)

assert(fe_above(fe_num(1), rmesh) == 401)
assert(fe_above(fe_num(401), rmesh) == 801)
assert(fe_above(fe_num(198*400+1), rmesh) == 199*400+1)
assert(fe_above(fe_num(200*400), rmesh) == 0)

assert(fe_below(fe_num(401), rmesh) == 1)
assert(fe_below(fe_num(400), rmesh) == 0)
assert(fe_below(fe_num(801), rmesh) == 401)
assert(fe_below(fe_num(199*400+1), rmesh) == 198*400+1)
assert(fe_below(fe_num(1), rmesh) == 0)


# Test finite element plane coordinates

# 10 rows x 20 cols mesh, with vertexes at each integer pair
mesh = RectMesh((0.,0.), (20.,10.), mesh_rc(10), mesh_rc(20))

assert(fe_coords(fe_num(1), mesh) == [0.0,0.0,1.0,1.0])
assert(fe_coords(fe_num(2), mesh) == [1.0,0.0,2.0,1.0])
assert(fe_coords(fe_num(20), mesh) == [19.0,0.0,20.0,1.0])
assert(fe_coords(fe_num(21), mesh) == [0.0,1.0,1.0,2.0])
assert(fe_coords(fe_num(200), mesh) == [19.0,9.0,20.0,10.0])
