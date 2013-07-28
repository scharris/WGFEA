using Base.Test

#using Winston

using WG
using Common
import Mesh, Mesh.fenum
import RMesh.RectMesh, RMesh.MeshCoord, RMesh.mesh_ldims, RMesh.mesh_ldim
import TMesh, TMesh.TriMesh
import WGBasis.WeakFunsPolyBasis
import VBF.AbstractVariationalBilinearForm
import VBF_a_s.a_s
import WGSol, WGSol.WGSolution
import Cubature.hcubature

function simple_2d_test_trimesh()
  u(x::Vector{R}) = x[1]*(1-x[1])*x[2]*(1-x[2])
  grad_u(x::Vector{R}) = [(1-2x[1])*x[2]*(1-x[2]), x[1]*(1-x[1])*(1-2x[2])]
  f(x::Vector{R}) = 2x[2]*(1-x[2]) + 2x[1]*(1-x[1])
  g = u

  ios = open("meshes/one_tri.msh")
  mesh = TriMesh(ios, 4)

  basis = WeakFunsPolyBasis(deg(3), deg(2), mesh)

  printBasisSummary(basis)

  println("Constructing VBF (computing vbf bel vs bel matrix)...")
  vbf = a_s(basis)
  println("VBF constructed.")
  wg = WGSolver(vbf, basis)

  println("Solving...")
  wg_sol = solve(f, g, wg)
  println("System solved.")

  # print_sample_points(wg_sol, u, grad_u)

  println("Computing L2 norm of error.")
  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(u, wg_sol))")

  #println("vbf semi-norm of Q u - u_h: $(WGSol.err_vs_proj_vbf_seminorm(u, wg_sol, vbf))")
  #  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(u, wg_sol))")
  #  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_grad_vs_wgrad_L2_norm(grad_u, wg_sol))")
end

function simple_2d_test()
  u(x::Vector{R}) = x[1]*(1-x[1])*x[2]*(1-x[2])
  grad_u(x::Vector{R}) = [(1-2x[1])*x[2]*(1-x[2]), x[1]*(1-x[1])*(1-2x[2])]
  f(x::Vector{R}) = 2x[2]*(1-x[2]) + 2x[1]*(1-x[1])
  g = 0.0

  mesh = RectMesh([0.,0.], [1.,1.], mesh_ldims(3,3))
  basis = WeakFunsPolyBasis(deg(2), deg(1), mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  println("Solving...")
  wg_sol = solve(f, g, wg)
  println("System solved.")

  # print_sample_points(wg_sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(u, wg_sol))")
  println("vbf semi-norm of Q u - u_h: $(WGSol.err_vs_proj_vbf_seminorm(u, wg_sol, vbf))")
  #  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(u, wg_sol))")
  #  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_grad_vs_wgrad_L2_norm(grad_u, wg_sol))")
end

function errs_and_diams_for_side_divs(errf::Function, space_dim::Dim, u::Function, f::Function, g::FunctionOrConst,
                                      mesh_bounds::(Vector{R}, Vector{R}), mesh_side_divs::Array{Int,1}, k::Deg)
  const num_runs = length(mesh_side_divs)

  const errs = Array(R, num_runs)
  const diams = Array(R, num_runs)
  for run=1:num_runs
    const mesh_ldims = fill(mesh_ldim(mesh_side_divs[run]+1), space_dim)
    const mesh = RectMesh(mesh_bounds[1], mesh_bounds[2], mesh_ldims)
    const basis = WeakFunsPolyBasis(k, deg(k-1), mesh)
    const vbf = a_s(basis)
    const wg = WGSolver(vbf, basis)
    const wg_sol = solve(f, g, wg)
    errs[run] = errf(u, wg_sol, vbf)
    diams[run] = Mesh.max_fe_diameter(mesh)
  end
  errs, diams
end


# Plot errors for approximation of function u where u(x) = cos(x_1) + sin(x_2), using the passed error function
# taking as arguments exact solution u, wg solution, and vbf.
function trig_2d_errs_plot(errf::Function, plot_name::String)
  u(x::Vector{R}) = cos(x[1]) + sin(x[2])
  # (grad u)(x) = (-sin(x_1), cos(x_2))
  # (div (grad u))(x) = -cos(x[1]) - sin(x[2])
  f(x::Vector{R}) = cos(x[1]) + sin(x[2])
  g(x::Vector{R}) = u(x)

  const mesh_bounds = ([0.,0.], [2pi,2pi])
  const integ_err_rel, integ_err_abs = 1e-12, 1e-12

  const mesh_side_divs = [20:5:35]

  const plot = FramedPlot()
  setattr(plot, "title", "Solving for u(x) = cos(x_1) + sin(x_2) on [0,2\\pi]^2")
  setattr(plot, "xlabel", "-log(h)")
  setattr(plot, "ylabel", "log ||Q_0 u - u_h||")

  const pt_types = ["diamond", "filled circle", "asterisk", "cross", "square", "triangle", "down-triangle", "right-triangle"]
  const pt_colors = ["blue", "green", "red", "cyan", "black", "yellow", "magenta"]

  const legend_pts_list = {}
  for k=deg(1):deg(4)
    const errs, diams = errs_and_diams_for_side_divs(errf, dim(2), u, f, g, mesh_bounds, mesh_side_divs, k)
    const num_errs = length(errs)
    const log_errs, minus_log_diams = Array(R, num_errs), Array(R, num_errs)
    for i=1:num_errs
      log_errs[i] = log(errs[i])
      minus_log_diams[i] = -log(diams[i])
    end

    const slope_str = @sprintf("%2.3f", -(log_errs[end]-log_errs[end-2])/(minus_log_diams[end]-minus_log_diams[end-2]))
    const pts = Points(minus_log_diams, log_errs, "type", pt_types[k], "color", pt_colors[k])
    setattr(pts, "label", "k=$(int(k)), O(h^{$slope_str})")
    const curve = Curve(minus_log_diams, log_errs, "color", pt_colors[k])
    add(plot, pts, curve)
    push!(legend_pts_list, pts)
  end

  const legend = Legend(.1,.25, legend_pts_list)
  add(plot, legend)

  file(plot, "plots/$plot_name.pdf")
end

function simple_4d_test()
  u(x::Vector{R}) = x[1]*(1-x[1]) * x[2]*(1-x[2]) * x[3]*(1-x[3]) * x[4]*(1-x[4])
  grad_u(x::Vector{R}) = [(1-2x[1]) * x[2]*(1-x[2]) * x[3]*(1-x[3]) * x[4]*(1-x[4]),
                          x[1]*(1-x[1]) * (1-2x[2]) * x[3]*(1-x[3]) * x[4]*(1-x[4]),
                          x[1]*(1-x[1]) * x[2]*(1-x[2]) * (1-2x[3]) * x[4]*(1-x[4]),
                          x[1]*(1-x[1]) * x[2]*(1-x[2]) * x[3]*(1-x[3]) * (1-2x[4])]
  f(x::Vector{R}) = 2*(x[2]*(1-x[2]) * x[3]*(1-x[3]) * x[4]*(1-x[4]) +
                       x[1]*(1-x[1]) * x[3]*(1-x[3]) * x[4]*(1-x[4]) +
                       x[1]*(1-x[1]) * x[2]*(1-x[2]) * x[4]*(1-x[4]) +
                       x[1]*(1-x[1]) * x[2]*(1-x[2]) * x[3]*(1-x[3]))
  g = 0.0

  mesh = RectMesh([0.,0.,0.,0.], [1.,1.,1.,1.], mesh_ldims(5,5,5,5))
  basis = WeakFunsPolyBasis(deg(3), deg(2), mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  wg_sol = solve(f, g, wg)

  print_sample_points(wg_sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(u, wg_sol))")
  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(u, wg_sol))")
  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_grad_vs_wgrad_L2_norm(grad_u, wg_sol))")
end


# u(x) = cos(|x|^2) in r^d
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
  f(x::Vector{R}) = let xx = dot(x,x); 2(d*sin(xx) + 2sum(i -> x[i]*x[i]*cos(xx), 1:d)) end
  g(x::Vector{R}) = u(x)

  mesh_min_bounds = zeros(R, d)
  mesh_max_bounds = ones(R, d)

  mesh = RectMesh(mesh_min_bounds, mesh_max_bounds, mesh_ldims)
  basis = WeakFunsPolyBasis(interior_polys_max_deg, side_polys_max_deg, mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  wg_sol = solve(f, g, wg)

  print_sample_points(wg_sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(u, wg_sol))")
  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(u, wg_sol))")
  println("vbf semi-norm of Q u - u_h: $(WGSol.err_vs_proj_vbf_seminorm(u, wg_sol, vbf))")
  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_grad_vs_wgrad_L2_norm(grad_u, wg_sol))")
end


function print_sample_points(wg_sol::WGSolution, u::Function, grad_u::Function)
  const d = Mesh.space_dim(basis.mesh)
  const int_rel_pt = zeros(R, d)
  const fe_origin = Array(R, d)

  for fe=fenum(1):Mesh.num_fes(basis.mesh)
    Mesh.fe_interior_origin!(fe, fe_origin, basis.mesh)
    fe_origin += int_rel_pt
    const wg_sol_val = WGSol.wg_sol_at_interior_rel_point(fe, int_rel_pt, wg_sol)
    const u_val = u(fe_or)

    const wgrad = WGSol.wg_sol_wgrad_on_interior(fe, wg_sol)
    const wgrad_val = Poly.value_at(wgrad, int_rel_pt)
    const grad_u_val = grad_u(fe_or)

    println("Point $fe_or: wg sol: $wg_sol_val, exact sol: $u_val")
    println("   wgrad: $wgrad_val, grad_u: $grad_u_val")
    println("   val diff: $(abs(wg_sol_val - u_val)),  grad diff: $(norm(wgrad_val - grad_u_val))")
  end
end


function printBasisSummary(basis::WeakFunsPolyBasis)
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  println("Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements.")
  println("Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println("Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")
  println("Data arrays for non-zeros are of initial (upper bound) size $(int64(interacting_bel_pairs_ub)).")
end

# Error functions
err_vs_proj_L2(u::Function, wg_sol::WGSolution, vbf::AbstractVariationalBilinearForm) =
  WGSol.err_vs_proj_L2_norm(u, wg_sol)

err_vs_proj_vbf_seminorm(u::Function, wg_sol::WGSolution, vbf::AbstractVariationalBilinearForm) =
  WGSol.err_vs_proj_vbf_seminorm(u, wg_sol, vbf)

#simple_2d_test()
simple_2d_test_trimesh()

#trig_Rd_test(mesh_ldims(5,5,5,5), deg(3), deg(2))

#simple_4d_test()

#trig_2d_errs_plot(err_vs_proj_vbf_seminorm, "err_vs_proj_vbf_seminorm")
#trig_2d_errs_plot(err_vs_proj_L2, "err_vs_proj_L2_norm")
