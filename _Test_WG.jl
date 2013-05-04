using Base.Test

using Winston

using WG
using Common
import Mesh, Mesh.fe_num
import RMesh.RectMesh, RMesh.MeshCoord, RMesh.mesh_ldims, RMesh.mesh_ldim
import WGBasis.WeakFunsPolyBasis
import VBF_a_s.a_s
import WGSol, WGSol.WGSolution
import Cubature.hcubature

# 75x75, deg(3),deg(2)
# non-cholesky solver: L2 norm of Q u - u_h: 4.694257769983945e-11 (repeated 4 times)
# cholesky solver: 1.169841805399121e-10 (repeated 3 times)

function simple_2d_test()
  u(x::Vector{R}) = x[1]*(1-x[1])*x[2]*(1-x[2])
  grad_u(x::Vector{R}) = [(1-2x[1])*x[2]*(1-x[2]), x[1]*(1-x[1])*(1-2x[2])]
  f(x::Vector{R}) = 2x[2]*(1-x[2]) + 2x[1]*(1-x[1])
  g = 0.0

  mesh = RectMesh([0.,0.], [1.,1.], mesh_ldims(75,75))
  basis = WeakFunsPolyBasis(deg(3), deg(2), mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  sol = solve(f, g, wg)

  # print_sample_points(sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(sol, u))")
  #  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(sol, u))")
  #  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_wgrad_vs_grad_L2_norm(sol, grad_u))")
  flush(OUTPUT_STREAM)
end

function errs_and_diams_for_side_divs(space_dim::Dim, u::Function, f::Function, g::FunctionOrConst,
                                      mesh_bounds::(Vector{R}, Vector{R}), mesh_side_divs::Array{Int,1}, k::Deg)
  const num_runs = length(mesh_side_divs)

  const l2_errs = Array(R, num_runs)
  const diams = Array(R, num_runs)
  for run=1:num_runs
    const mesh_ldims = fill(mesh_ldim(mesh_side_divs[run]+1), space_dim)
    const mesh = RectMesh(mesh_bounds[1], mesh_bounds[2], mesh_ldims)
    const basis = WeakFunsPolyBasis(k, deg(k-1), mesh)
    const vbf = a_s(basis)
    const wg = WGSolver(vbf, basis)
    const sol = solve(f, g, wg)
    l2_errs[run] = WGSol.err_vs_proj_L2_norm(sol, u)
    diams[run] = Mesh.max_fe_diameter(mesh)
  end
  l2_errs, diams
end


# Plot errors for approximation of function u where u(x) = cos(x_1) + sin(x_2).
function trig_2d_errs_plot()
  u(x::Vector{R}) = cos(x[1]) + sin(x[2])
  # (grad u)(x) = (-sin(x_1), cos(x_2))
  # (div (grad u))(x) = -cos(x[1]) - sin(x[2])
  f(x::Vector{R}) = cos(x[1]) + sin(x[2])
  g(x::Vector{R}) = u(x)

  const mesh_bounds = ([0.,0.], [4pi,4pi])
  const integ_err_rel, integ_err_abs = 1e-14, 1e-14

  const mesh_side_divs = [20:5:65]

  const plot = FramedPlot()
  setattr(plot, "title", "Solving for u(x) = cos(x_1) + sin(x_2) on [0,4\\pi]^2")
  setattr(plot, "xlabel", "-log(h)")
  setattr(plot, "ylabel", "log ||Q_0 u - u_h||")

  const pt_types = ["diamond", "filled circle", "asterisk", "cross", "square", "triangle", "down-triangle", "right-triangle"]
  const pt_colors = ["blue", "green", "red", "cyan", "black", "yellow", "magenta"]

  const legend_pts_list = {}
  for k=deg(1):deg(4)
    const errs, diams = errs_and_diams_for_side_divs(dim(2), u, f, g, mesh_bounds, mesh_side_divs, k)
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

  file(plot, "plots/errs_plot.pdf")
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

  sol = solve(f, g, wg)

  print_sample_points(sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(sol, u))")
  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(sol, u))")
  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_wgrad_vs_grad_L2_norm(sol, grad_u))")
  flush(OUTPUT_STREAM)
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
  f(x::Vector{R}) = let xx = dot(x,x) 2(d*sin(xx) + 2sum(i -> x[i]*x[i]*cos(xx), 1:d)) end
  g(x::Vector{R}) = u(x)

  mesh_min_bounds = zeros(R, d)
  mesh_max_bounds = ones(R, d)

  mesh = RectMesh(mesh_min_bounds, mesh_max_bounds, mesh_ldims)
  basis = WeakFunsPolyBasis(interior_polys_max_deg, side_polys_max_deg, mesh)
  vbf = a_s(basis)
  wg = WGSolver(vbf, basis)

  wg_sol = solve(f, g, wg)

  print_sample_points(wg_sol, u, grad_u)

  println("L2 norm of Q u - u_h: $(WGSol.err_vs_proj_L2_norm(wg_sol, u))")
  println("L2 norm of u - u_h: $(WGSol.err_L2_norm(wg_sol, u))")
  println("L2 norm of grad u - wgrad u_h: $(WGSol.err_wgrad_vs_grad_L2_norm(wg_sol, grad_u))")
  flush(STDOUT)
end


function print_sample_points(wg_sol::WGSolution, u::Function, grad_u::Function)
  const int_rel_pt = zeros(R, Mesh.space_dim(basis.mesh))
  for fe=fe_num(1):Mesh.num_fes(basis.mesh)
    const fe_or = Mesh.fe_interior_origin(fe, basis.mesh) + int_rel_pt
    const wg_sol_val = WGSol.wg_sol_at_interior_rel_point(wg_sol, fe, int_rel_pt)
    const u_val = u(fe_or)

    const wgrad = WGSol.wg_sol_wgrad_on_interior(wg_sol, fe)
    const wgrad_val = Poly.value_at(wgrad, int_rel_pt)
    const grad_u_val = grad_u(fe_or)

    println("Point $fe_or: wg sol: $wg_sol_val, exact sol: $u_val")
    println("   wgrad: $wgrad_val, grad_u: $grad_u_val")
    println("   val diff: $(abs(wg_sol_val - u_val)),  grad diff: $(norm(wgrad_val - grad_u_val))")
  end
  flush(OUTPUT_STREAM)
end

simple_2d_test()

#trig_Rd_test(mesh_ldims(5,5,5,5), deg(3), deg(2))

#simple_4d_test()

# trig_2d_errs_plot()

trig_2d_sob_norms_unused = "
  # (||u||_2)^2 = int_{mesh_bounds} u(x)^2 + |(grad u)(x)|^2  + (D_1,1 u)(x)^2 + (D_2,2 u)(x)^2 + (D_1,2 u)(x)^2 dx
  sob2_norm_sq_u = begin
    function sum_sq_partials(x::Vector{R})
      const cos_x1, cos_x2, sin_x1, sin_x2 = cos(x[1]), cos(x[2]), sin(x[1]), sin(x[2])
      const u_sq = let s = cos_x1 + sin_x2  s*s end
      # (D_1 u)(x) = -sin(x_1),   (D_2 u)(x) = cos(x_2)
      const grad_u_norm_sq = sin_x1 * sin_x1 + cos_x2 * cos_x2
      const d11_sq = cos_x1 * cos_x1 # (D_1,1 u)(x) = -cos(x_1)
      const d22_sq = sin_x2 * sin_x2 # (D_2,2 u)(x) = -sin(x_2)
      const d12_sq = 0.
      u_sq + grad_u_norm_sq + (d11_sq + d22_sq + d12_sq)
    end
    hcubature(sum_sq_partials, mesh_bounds[1], mesh_bounds[2], integ_err_rel, integ_err_abs)[1]
  end

  sob3_norm_sq_u = begin
    function sum_additional_sq_partials(x::Vector{R})
      const cos_x2, sin_x1 = cos(x[2]), sin(x[1])
      const d111_sq = sin_x1 * sin_x1 # (D_1,1,1 u)(x) = sin(x_1)
      const d222_sq = cos_x2 * cos_x2 # (D_2,2,2 u)(x) = -cos(x_2)
      d111_sq + d222_sq # Mixed partials are 0.
    end
    sob2_norm_sq_u + hcubature(sum_additional_sq_partials, mesh_bounds[1], mesh_bounds[2], integ_err_rel, integ_err_abs)[1]
  end

  sob4_norm_sq_u = begin
    function sum_additional_sq_partials(x::Vector{R})
      const cos_x1, sin_x2 = cos(x[1]), sin(x[2])
      const d1111_sq = cos_x1 * cos_x1 # (D_1,1,1,1 u)(x) = cos(x_1)
      const d2222_sq = sin_x2 * sin_x2 # (D_2,2,2,2 u)(x) = sin(x_2)
      d1111_sq + d2222_sq # Mixed partials are 0.
    end
    sob3_norm_sq_u + hcubature(sum_additional_sq_partials, mesh_bounds[1], mesh_bounds[2], integ_err_rel, integ_err_abs)[1]
  end

  sob5_norm_sq_u = begin
    function sum_additional_sq_partials(x::Vector{R})
      const sin_x1, cos_x2 = sin(x[1]), cos(x[2])
      const d11111_sq = sin_x1 * sin_x1 # (D_1,1,1,1,1 u)(x) = -sin(x_1)
      const d22222_sq = cos_x2 * cos_x2 # (D_2,2,2,2,2 u)(x) = cos(x_2)
      # All mixed partials are 0 as before.
      d11111_sq + d22222_sq
    end
    sob4_norm_sq_u + hcubature(sum_additional_sq_partials, mesh_bounds[1], mesh_bounds[2], integ_err_rel, integ_err_abs)[1]
  end

  const u_sob_norms = [NaN, sqrt(sob2_norm_sq_u), sqrt(sob3_norm_sq_u), sqrt(sob4_norm_sq_u), sqrt(sob5_norm_sq_u)]
"
