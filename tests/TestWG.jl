using Test
require("../Common")
require("../Poly")
require("../Mesh")
require("../WGrad")
require("../WGBasis")
require("../RMesh")
require("../Proj")
require("../VBF")
require("../WG")


using WG

using Common
import RMesh.RectMesh, RMesh.mesh_coord

mesh = RectMesh([0.,0.], [3.,3.], [mesh_coord(3), mesh_coord(3)])
bf = x -> x

wg_solver = WGSolver(bf, deg(4), mesh)
