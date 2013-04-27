using Base.Test
using WG
using Common
import Mesh, Mesh.fe_num
import RMesh.RectMesh, RMesh.MeshCoord, RMesh.mesh_dims
import WGBasis.WeakFunsPolyBasis
import VBF_a_s.a_s
import VBF_a.a
import Sol, Sol.WGSolution


function simple_2d_test()
  u(x::Vector{R}) = x[1]*(1-x[1])*x[2]*(1-x[2])
  grad_u(x::Vector{R}) = [(1-2x[1])*x[2]*(1-x[2]), x[1]*(1-x[1])*(1-2x[2])]
  f(x::Vector{R}) = 2x[2]*(1-x[2]) + 2x[1]*(1-x[1])
  g = 0.0

  mesh = RectMesh([0.,0.], [1.,1.], mesh_dims(5,5))
  basis = WeakFunsPolyBasis(deg(2), deg(1), mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  sol = solve(f, g, wg)

  #print_sample_points(sol, u, grad_u, basis)

  println("L2 norm of u - u_h: $(Sol.err_L2_norm(sol, u, basis))")
  println("L2 norm of grad u - wgrad u_h: $(Sol.err_wgrad_vs_grad_L2_norm(sol, grad_u, basis))")
  flush(OUTPUT_STREAM)
end

# u(x) = cos(|x|^2) in R^d
function trig_Rd_test(mesh_ldims::Array{MeshCoord,1},
                      interior_polys_max_deg::Deg,
                      side_polys_max_deg::Deg)
  const d = length(mesh_ldims)

  u(x::Vector{R}) = cos(dot(x,x))
  # ((grad u)(x))_i = -2 sin(|x|^2) x_i
  grad_u(x::Vector{R}) = -2sin(dot(x,x))*x
  # (D_i (grad u)_i)(x) = -2 ( cos(|x|^2)(2 x_i) x_i + sin(|x|^2) )
  # So -(div grad u)(x) =  2 sum_{i=1..d}{ 2 (x_i)^2 cos(|x|^2)  + sin(|x|^2) }
  #                     =  2 ( d sin(|x|^2) + 2 sum_{i=1..d}{ (x_i)^2 cos(|x|^2) } )
  f(x::Vector{R}) = let xx = dot(x,x) 2(d*sin(xx) + 2sum(i -> x[i]*x[i]*cos(xx), 1:d)) end
  g(x::Vector{R}) = u(x)

  mesh_min_bounds = zeros(R, d)
  mesh_max_bounds = ones(R, d)

  mesh = RectMesh(mesh_min_bounds, mesh_max_bounds, mesh_ldims)
  basis = WeakFunsPolyBasis(interior_polys_max_deg, side_polys_max_deg, mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  println(STDERR, "WGSolver created, $(int64(wg.basis.total_bels)) basis elements, $(nnz(wg.vbf_bel_vs_bel_transpose)) nonzeros in vbf bel vs bel matrix.")
  wg_sol = solve(f, g, wg)
  println(STDERR, "Solve completed.")

  print_sample_points(wg_sol, u, grad_u, basis)

  println(STDERR, "L2 norm of u - u_h: $(Sol.err_L2_norm(wg_sol, u, basis))")
  println(STDERR, "L2 norm of grad u - wgrad u_h: $(Sol.err_wgrad_vs_grad_L2_norm(wg_sol, grad_u, basis))")
  flush(OUTPUT_STREAM)
end


function print_sample_points(wg_sol::WGSolution, u::Function, grad_u::Function, basis::WeakFunsPolyBasis)
  for fe=fe_num(1):Mesh.num_fes(basis.mesh)
    const fe_or = Mesh.fe_interior_origin(fe, basis.mesh)
    const wg_sol_val = Sol.wg_sol_at_interior_rel_point(wg_sol, fe, [.0,.0], basis)
    const u_val = u(fe_or)

    const wgrad = Sol.wg_sol_wgrad_on_interior(wg_sol, fe, basis)
    const wgrad_val = Poly.value_at(wgrad, [.0,.0])
    const grad_u_val = grad_u(fe_or)

    println(STDERR, "wg sol: $wg_sol_val, exact sol: $u_val, wgrad: $wgrad_val, grad_u: $grad_u_val")
    println(STDERR, "  ==> val diff: $(abs(wg_sol_val - u_val)),  grad diff: $(norm(wgrad_val - grad_u_val))")
  end
  flush(OUTPUT_STREAM)
end

#simple_2d_test()

trig_Rd_test(mesh_dims(5,5,5,5), deg(3), deg(2))
# "[Finished in 2809.1s]"
