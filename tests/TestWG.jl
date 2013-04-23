using Test
require("../Common")
require("../Poly")
require("../Mesh")
require("../RMesh")
require("../WGrad")
require("../WGBasis")
require("../Proj")
require("../VBF")
require("../VBF_a_s")
require("../Sol")
require("../WG")

using WG
using Common
import Mesh, Mesh.fe_num
import RMesh.RectMesh, RMesh.mesh_coord
import WGBasis.WeakFunsPolyBasis
import VBF_a_s.a_s
import Sol

f(x::Vector{R}) = 2x[1]*(1-x[1]) + 2x[2]*(1-x[2])
g(x::Vector{R}) = 0

mesh = RectMesh([0.,0.], [1.,1.], [mesh_coord(5), mesh_coord(5)])
basis = WeakFunsPolyBasis(deg(4), deg(3), mesh)
vbf = a_s(basis)
wg = WGSolver(vbf, basis)

sol = solve(f, g, wg)

u(x::Vector{R}) = x[1]*(1-x[1])*x[2]*(1-x[2])
grad_u(x::Vector{R}) = [(1-2x[1])*x[2]*(1-x[2]), x[1]*(1-x[1])*(1-2x[2])]

for fe=fe_num(1):Mesh.num_fes(mesh)
  const fe_or = Mesh.fe_interior_origin(fe, mesh)
  const sol_val = Sol.solution_value_at_interior_rel_point(sol, fe, [.0,.0], basis)
  const u_val = u(fe_or)

  const wgrad = Sol.solution_wgrad_on_interior(sol, fe, basis)
  const wgrad_val = Poly.value_at(wgrad, [.0,.0])
  const grad_u_val = grad_u(fe_or)

  println("sol: $sol_val, actual: $u_val, wgrad: $wgrad_val, grad_u: $grad_u_val")
  println("  ==> val diff: $(abs(sol_val - u_val)),  grad diff: $(norm(wgrad_val - grad_u_val))")
end

err_L2 = Sol.err_est_L2(sol, u, basis)
println("-----------------------------------------------")
println("L2 norm of difference: $(err_L2)")
println("-----------------------------------------------")
norm(a) = sqrt(dot(a,a))
