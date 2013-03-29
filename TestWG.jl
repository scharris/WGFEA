using Test
using WG

using Common
import RMesh2.RectMesh2

mesh = RectMesh2((0.,0.), (3.,3.), 3, 3)
bf = x -> x

wg_solver = WGSolver(bf, deg(4), mesh)
