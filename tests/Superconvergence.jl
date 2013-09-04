
require("Common.jl")
require("Poly.jl")
require("Mesh.jl")
require("RMesh.jl")
require("TMesh.jl")
require("WGBasis.jl")
require("VBF.jl")
require("VBF_a_s.jl")
require("WG.jl")
require("WGSol.jl")

using Common
import Poly, Poly.Monomial, Poly.Polynomial, Poly.grad
import Mesh, Mesh.AbstractMesh, Mesh.FENum, Mesh.fenum
import RMesh, RMesh.RectMesh, RMesh.MeshCoord, RMesh.mesh_coord, RMesh.fe_mesh_coords!
import TMesh, TMesh.TriMesh, TMesh.Point, TMesh.BaseTri
import WGBasis, WGBasis.WeakFunsPolyBasis
import VBF_a_s.a_s
import WG, WG.WGSolver
import WGSol, WGSol.WGSolution
import Cubature.hcubature

# Method
# The method is to let the coarser, tau element diameter mesh be a rectangular mesh formed by
# dividing a rectangular region into a grid of equal rectangles, and to have the finer h element
# diameter mesh consist of triangles formed by subdividing the rectangles of the tau mesh by some
# number of iterations.  Thus every tau mesh element is the union of some subset of h mesh
# elements, which makes projection tractable with reasonable efficiency.
#
# Let Omega be the rectangular region to be meshed, and let alpha in (0,1) be rational.
# Our meshing strategy requires that there are positive integers n and m such that:
#   1) tau/h = 2^n,
#      because the triangles of diameter h are formed by subdividing the tau mesh rectangles, with
#      initial hypotenuses along the rectangle diagonals prior to subdivision, and
#   2) tau = diam(Omega)/m,
#      because the tau mesh is formed by evenly dividing the original rectangular regionsin into an
#      m x m grid, with resulting rectangles having having diameter diam(Omega)/m.
#   3) tau = h^alpha
# These conditions together imply that n and m satisfy
#  (4)  2^{n alpha/(alpha-1)} = diam(Omega)/m.
# Conversely, if postive integers m and n satisfy (4), then
# by letting tau_{m,n} = 2^{n alpha/(alpha-1)} and h_{m,n} = tau_{m,n}^{1/alpha},
# then (1) and (3) follow directly, and (2) follows from
#   tau_{m,n}/h = tau_{m,n}^{1 - 1/alpha} = tau_{m,n}^{(alpha-1)/alpha} = 2^n .
# Thus (4) implies (1), (2), and (3).
#
# We now only need to provide a sequence of pairs {(n_i,m_i): i = 1,2,...} satisfying (4), where the n_i are
# increasing so that h->0 as i->infinity.
#
# To do this, assume that we've found a single pair n and m satisyfing (4).
# Since alpha is rational, there are integers p < 0 and q > 0 such that 
#    alpha/(alpha-1) = p/q  .
# Then
#    2^{(n+q) alpha/(alpha-1)} = 2^{n alpha/(alpha-1)} 2^p
#                              = (diam(Omega)/m) 2^p 
#                              = diam(Omega)/(m 2^|p|) (since p < 0)
#                              
# Thus we see that if m,n satisfy our condition (4), then so do the integers m 2^|p|, n+q in place of m
# and n.   This gives us our sequence of pairs,
#   {(m 2^{j p}, n + j q), j=0,1,...},
# of mesh parameters satisfying (1),(2), and (3) and with h->0, so long as our initial pair n,m satisfies
# considion (4).
#
# Example
# -------
# Let s = 2, k = 2, r = 4.
# Then alpha = 3/5 and alpha/(alpha-1) = -3/2.
# By (4), we need to find n and m (and diam(Omega)) with:
#    2^{n alpha/(alpha-1)} = diam(Omega)/m
#
# So taking log2 of both sides, in this case we need
#   (-3/2) n = log2(diam(Omega)) - log2(m)
#   n = (-2/3) (log2(diam(Omega)) - log2(m))
# Choose the region Omega so that diam(Omega) = 1
# Then we need
#   n = (2/3) log2(m)
# We can let m = 8, and n = 2 to satisfy this.
#
# Then to verify:
# tau is the region diameter over an integer:
#    tau = 1/8 = diam(Omega)/8,
# and setting h = tau^(1/alpha) = 0.03125,
# we have
#   tau/h = 4 = 2^2
# so the h mesh can be obtained by subdivision of the tau mesh.
# And of course tau = h^alpha by choice of h.

# The next mesh parameters in the sequence after (m,n) are (2^|p| m, n + q) = (8 m, n + 2) = (64, 4).
# Then tau = 1/64, h = tau^(1/alpha) = 0.0009765625, and 
# tau/h = 16 = 2^4.

immutable LaplaceModelProblem
  laplace_eq_rhs::Function
  boundary_fn::Function
  region_min::Array{R,1}
  region_max::Array{R,1}
  sol_H_regularity::Rational{Int}
  sol::Function
  grad_sol::Array{Function,1}
end

immutable MeshPairSpec
  tau_mesh_divs_per_side::Int # m
  h_mesh_subdiv_ops::Int      # n
end

immutable ErrorEst
  max_fe_diam::R
  err_wg_L2::R
  err_wg_H1::R
  err_proj_L2::R
  err_proj_H1::R
end


const integration_rel_err = 1e-12
const integration_abs_err = 1e-12

function test_superconvergence(mod_prob::LaplaceModelProblem,
                               wg_polys_max_deg::Deg,
                               proj_space_polys_max_deg::Deg,
                               mesh_pair_specs_list::Array{MeshPairSpec})
  const s = mod_prob.sol_H_regularity
  const k = wg_polys_max_deg
  const r = proj_space_polys_max_deg
  const alpha = (k + s - 1)//(r + 1 - min(0, 2 - s))

  function solve_and_project(mesh_pair_spec)
    const m = mesh_pair_spec.tau_mesh_divs_per_side
    const tau = norm(mod_prob.region_max - mod_prob.region_min) / m
    const tau_mesh = RectMesh(mod_prob.region_min, mod_prob.region_max, [mesh_coord(m), mesh_coord(m)])
    
    const n = mesh_pair_spec.h_mesh_subdiv_ops
    const h_mesh = TMesh.from_rect_mesh(tau_mesh, n)
    const h = tau^(1/alpha)
    
    println(h_mesh)
    printMeshSizesInfo(h_mesh, tau_mesh, tau, alpha, n)

    const wg_sol = wg_solution(mod_prob, h_mesh, k)

    const projs_by_tau_mesh_fe = project(wg_sol, tau_mesh, r)

    const err_wg_L2 = WGSol.err_L2_norm(mod_prob.sol, wg_sol)
    const err_wg_H1 = WGSol.err_H1_norm(mod_prob.sol, mod_prob.grad_sol, wg_sol)
    const err_proj_L2 = err_proj_L2_norm(mod_prob.sol, projs_by_tau_mesh_fe, tau_mesh)
    const err_proj_H1 = err_proj_H1_norm(mod_prob.sol, mod_prob.grad_sol, projs_by_tau_mesh_fe, tau_mesh)

    ErrorEst(h, err_wg_L2, err_wg_H1, err_proj_L2, err_proj_H1)
  end

  const map_fn = ParCtrl.parallel_superconvergence_solve_and_project ? pmap : map
  const errs = map_fn(solve_and_project, mesh_pair_specs_list)
  
  println("")
  println("Convergence Results for s=$s, k=$k, r=$r")
  println("----------------------------------------------")
  println("h, err_wg_L2, err_wg_H1, err_proj_L2, err_proj_H1")

  for i=1:length(errs)
    const err = errs[i]
    println("$(err.max_fe_diam), $(err.err_wg_L2), $(err.err_wg_H1), $(err.err_proj_L2), $(err.err_proj_H1)")
  end

  println("----------------------------------------------")
  println("h1, h2, wg_rate_L2, wg_rate_H1, proj_rate_L2, proj_rate_H1")
  
  for i=1:length(errs)-1
    const err1 = errs[i]
    const err2 = errs[i+1]
    const slope_wg_L2 = (log(err2.err_wg_L2)-log(err1.err_wg_L2)) / (log(err2.max_fe_diam)-log(err1.max_fe_diam))
    const slope_wg_H1 = (log(err2.err_wg_H1)-log(err1.err_wg_H1)) / (log(err2.max_fe_diam)-log(err1.max_fe_diam))
    const slope_proj_L2 = (log(err2.err_proj_L2)-log(err1.err_proj_L2)) / (log(err2.max_fe_diam)-log(err1.max_fe_diam))
    const slope_proj_H1 = (log(err2.err_proj_H1)-log(err1.err_proj_H1)) / (log(err2.max_fe_diam)-log(err1.max_fe_diam))
    println("$(err1.max_fe_diam), $(err2.max_fe_diam), $slope_wg_L2, $slope_wg_H1, $slope_proj_L2, $slope_proj_H1")
  end
end

function wg_solution(mod_prob::LaplaceModelProblem, mesh::AbstractMesh, k::Deg)
  println("Constructing basis.")
  const basis = WGBasis.WeakFunsPolyBasis(k, deg(k-1), mesh)
  printBasisSummary(basis)

  println("Constructing VBF.")
  const vbf = a_s(basis)

  println("Constructing WG solver and computing vbf bel vs bel matrix.")
  const wg = WG.WGSolver(vbf, basis)

  println("Solving...")
  const wg_sol = WG.solve(mod_prob.laplace_eq_rhs, mod_prob.boundary_fn, wg)
  println("System solved.")

  wg_sol
end


#  Project the wg sol onto the space of piecewise polynomials of deg <= r on the passed rectangular mesh.
#  Returns an array of polynomials, indexed by tau mesh element number.
#  Method
#  Given an element e_tau of the tau mesh, and a basis {v_i} of P_r(e_tau), we have
#    (Q_tau u_h, v_i)_{e_tau} = (u_h, v_i)_{e_tau} for i=1..dim(P_r(e_tau)).
#  Letting e_{h,j} be the h-mesh elements whose union is e_tau (which exist by construction),
#  this gives us
#    (Q_tau u_h, v_i)_{e_tau} = sum_j (u_h, v_i)_{e_{h,j}}.
#  Letting lambda be the vector of unknown coefficients satisfying
#    Q_tau u_h = sum_j lambda_j v_j,
#  we obtain the linear system
#  (sys)  M lambda = b
#  where M_ij = (v_j, v_i)_{e_tau},
#    and b_i = sum_j (u_h, v_i)_{e_{h,j}} for i=1..dim(P_r(e_tau)).
function project(wg_sol::WGSolution, tau_mesh::RectMesh, r::Deg)
  println("Projecting wg solution.")
  const h_mesh = wg_sol.basis.mesh
  const num_tau_mesh_fes = Mesh.num_fes(tau_mesh)
  const included_h_fes_by_tau_fe = fenums_by_geoment(h_mesh, num_tau_mesh_fes)

  # Reference basis monomials, which when interpreted locally on a tau mesh element will form
  # a basis of the space onto which we are projecting for the element.
  const proj_space_basis_mons = Poly.mons_of_deg_le(r, Mesh.space_dim(tau_mesh))                   # v_*
  const proj_space_dim = length(proj_space_basis_mons)
  const proj_space_basis_vs_basis_ips = ips_matrix(proj_space_basis_mons, RMesh.fe_dims(tau_mesh)) # M

  function project_onto(tau_fe::FENum)
    const sys_rhs = zeros(R, proj_space_dim)                                                       # b
    # Compute rhs (b) column vector of (sys).
    for i=1:proj_space_dim # (sys) equation number, index into proj space basis
      # Compute component i of (sys) rhs.
      const tau_mon = proj_space_basis_mons[i]                                                     # v_i
      const h_fes = included_h_fes_by_tau_fe[tau_fe]                                               # e_{h,*}
      for h_fe in h_fes                                                                            # e_{h,j}
        const u_h_rel = WGSol.wg_sol_interior_poly(h_fe, wg_sol)
        # Translate the tau monomial to a polynomial to be interpreted relative to the h_fe origin.
        const origin_diff = Mesh.fe_interior_origin(h_fe, h_mesh) - Mesh.fe_interior_origin(tau_fe, tau_mesh)
        const tau_mon_h_fe_rel = Poly.canonical_form(Poly.translate(tau_mon, origin_diff))
        const u_h_x_tau_mon = Poly.canonical_form(u_h_rel * tau_mon_h_fe_rel)
        const oshape = Mesh.oriented_shape_for_fe(h_fe, h_mesh)
        sys_rhs[i] += Mesh.integral_face_rel_on_oshape_face(u_h_x_tau_mon, oshape, Mesh.interior_face, h_mesh)
      end
    end
    const proj_coefs = proj_space_basis_vs_basis_ips \ sys_rhs
    Polynomial(proj_space_basis_mons, proj_coefs)
  end

  const map_fn = ParCtrl.parallel_superconvergence_project ? pmap : map
  const projs = map_fn(project_onto, [fenum(1):num_tau_mesh_fes])

  convert(Array{Polynomial,1}, projs)
end


function fenums_by_geoment(tri_mesh::TriMesh, num_geoments::FENum)
  const fenums_by_ge = sizehint(Dict{TMesh.Tag, Array{FENum,1}}(), num_geoments)
  const approx_fes_per_geoment = int64(Mesh.num_fes(tri_mesh) / num_geoments)
  for fe=fenum(1):Mesh.num_fes(tri_mesh)
    const ge = TMesh.geometric_entity_tag(fe, tri_mesh)
    ge_fes = get(fenums_by_ge, ge, nothing)
    if ge_fes == nothing
      ge_fes = sizehint(Array(FENum, 0), int(approx_fes_per_geoment * 1.10))
      fenums_by_ge[ge] = ge_fes
    end
    push!(ge_fes, fe)
  end
  fenums_by_ge
end


function ips_matrix(mons::Array{Monomial,1}, dom_rect_dims::Array{R,1})
  const num_mons = length(mons)
  const m = Array(R, num_mons,num_mons)
  for i=1:num_mons
    const mon_i = mons[i]
    for j=1:i-1
      const ip = Poly.integral_on_rect_at_origin(mon_i * mons[j], dom_rect_dims)
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = Poly.integral_on_rect_at_origin(mon_i * mon_i, dom_rect_dims)
  end
  m
end


function err_proj_L2_norm(exact_sol::Function,
                          approx_sol_by_rmesh_fe::Array{Polynomial,1},
                          rmesh::RectMesh)
  const d = Mesh.space_dim(rmesh)
  const zs = zeros(R,d)
  const fe_dims = RMesh.fe_dims(rmesh)

  # Work data for sq_err function below.
  const x = Array(R, d)
  const fe_origin = Array(R, d)

  integral_sq_err = zeroR
  for fe=fenum(1):Mesh.num_fes(rmesh)
    Mesh.fe_interior_origin!(fe, fe_origin, rmesh)
    
    function sq_err(fe_rel_x::Vector{R})
      for i=1:d
        x[i] = fe_rel_x[i] + fe_origin[i]
      end
      const err = exact_sol(x) - Poly.value_at(approx_sol_by_rmesh_fe[fe], fe_rel_x)
      err * err
    end
    
    integral_sq_err += hcubature(sq_err,
                                 zs, fe_dims,
                                 reltol=integration_rel_err, abstol=integration_abs_err)[1]
  end
  sqrt(integral_sq_err)
end


function err_proj_H1_norm(exact_sol::Function,
                          grad_exact_sol::Array{Function,1},
                          approx_sol_by_rmesh_fe::Array{Polynomial,1},
                          rmesh::RectMesh)
  const d = Mesh.space_dim(rmesh)
  const zs = zeros(R,d)
  const fe_dims = RMesh.fe_dims(rmesh)

  const grad_approx_sol_by_rmesh_fe = map(grad, approx_sol_by_rmesh_fe)

  # Work data for sum_sq_errs function below.
  const x = Array(R, d)
  const fe_origin = Array(R, d)

  integral_sq_errs = zeroR

  for fe=fenum(1):Mesh.num_fes(rmesh)
    Mesh.fe_interior_origin!(fe, fe_origin, rmesh)
    const approx = approx_sol_by_rmesh_fe[fe]
    const grad_approx = grad_approx_sol_by_rmesh_fe[fe]
    
    function sum_sq_errs(fe_rel_x::Vector{R})
      # set the global point coordinates
      for i=1:d
        x[i] = fe_rel_x[i] + fe_origin[i]
      end
      # compute the sum of the square first partial derivative errors
      sum_first_partial_sq_errs = zeroR
      for n=dim(1):d
        const partial_diff = grad_exact_sol[n](x) - Poly.value_at(grad_approx[n], fe_rel_x)
        sum_first_partial_sq_errs += partial_diff * partial_diff
      end
      const err = exact_sol(x) - Poly.value_at(approx, fe_rel_x)
      err * err + sum_first_partial_sq_errs
    end
    
    integral_sq_errs += hcubature(sum_sq_errs,
                                  zs, fe_dims,
                                  reltol=integration_rel_err, abstol=integration_abs_err)[1]
  end

  sqrt(integral_sq_errs)
end

function printBasisSummary(basis::WeakFunsPolyBasis)
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  println("Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements.")
  println("Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println("Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")
  println("Data arrays for non-zeros are of initial (upper bound) size $(int64(interacting_bel_pairs_ub)).")
end

function printMeshSizesInfo(h_mesh, tau_mesh, tau, alpha, tri_subdivs)
  const h = tau^(1/alpha) 
  const h_act = Mesh.max_fe_diameter(h_mesh)
  const tau_act = Mesh.max_fe_diameter(tau_mesh)
  println()
  println("h=$h_act, tau=$tau_act, tau/h^alpha=$(tau_act/h_act^alpha), tau/h=$(tau_act/h_act)")
end



function makeMeshPairSpecs(init_mesh_pair_spec::MeshPairSpec, numMeshes::Int, alpha::Rational{Int})
  const meshSpecs = Array(MeshPairSpec, numMeshes)
  const m = init_mesh_pair_spec.tau_mesh_divs_per_side
  const n = init_mesh_pair_spec.h_mesh_subdiv_ops
  const alpha_over_alpha_minus_1 = alpha/(alpha-1)
  const p = abs(alpha_over_alpha_minus_1.num)
  const q = abs(alpha_over_alpha_minus_1.den)
  for j=0:numMeshes-1
    meshSpecs[j+1] = MeshPairSpec(m * 2^(j*p), n + j*q)
  end
  meshSpecs
end

# function go_orig()
#   u(x::Vector{R}) = cos(x[1]) + sin(x[2])
#   grad_u = [x -> -sin(x[1]), x -> cos(x[2])]
#   # (div (grad u))(x) = -cos(x[1]) - sin(x[2])
#   f(x::Vector{R}) = cos(x[1]) + sin(x[2])
#   g(x::Vector{R}) = u(x)
#   const s = 2//1
#   const mod_prob = LaplaceModelProblem(f, g, [0.,0.], [1./sqrt(2), 1./sqrt(2)], s, u, grad_u)
#   const k = deg(1)
#   const r = deg(3)
#   const alpha = (k + s - 1)/(r + 1 - min(0, 2 - s))

#   println("Expected L2 convergence rate is $((k + s - 1)*(r + 1)/(r + 1 + max(0, s-2))).")
#   println("Expected H1 convergence rate is $((k + s - 1)*r/(r + 1 + max(0, s-2))).")

#   const num_mesh_pairs = 3
#   const mesh_pair_specs = makeMeshPairSpecs(MeshPairSpec(2,1), num_mesh_pairs, alpha)
#   test_superconvergence(mod_prob,
#                         k,
#                         r,
#                         mesh_pair_specs)
# end

function go()
  ############## example 1
  #u(x::Vector{R}) = cos(x[1]) + sin(x[2])
  #grad_u = [x -> -sin(x[1]), x -> cos(x[2])]
  # (div (grad u))(x) = -cos(x[1]) - sin(x[2])
  #f(x::Vector{R}) = cos(x[1]) + sin(x[2])
  #g(x::Vector{R}) = u(x)
  ############## example 2
  #u(x::Vector{R}) = x[1]*x[2]*(x[1]-1./sqrt(2))*(x[2]-1./sqrt(2))
  #grad_u = [x -> x[2]*(x[2]-1/sqrt(2))*(2.*x[1]-1./sqrt(2)), x -> x[1]*(x[1]-1/sqrt(2))*(2.*x[2]-1./sqrt(2))]
  # (div (grad u))(x) = 2.*x[2]*(x[2]-1/sqrt(2))+2.*x[1]*(x[1]-1/sqrt(2))
  #f(x::Vector{R}) = -2.*x[2]*(x[2]-1/sqrt(2))-2.*x[1]*(x[1]-1/sqrt(2))
  #g(x::Vector{R}) = u(x)
  ############## example 3
  #u(x::Vector{R}) = sin(sqrt(2)*pi*x[1])*cos(sqrt(2)*pi*x[2])
  #grad_u = [x -> sqrt(2)*pi*cos(sqrt(2)*pi*x[1])*cos(sqrt(2)*pi*x[2]), x -> -sqrt(2)*pi*sin(sqrt(2)*pi*x[1])*sin(sqrt(2)*pi*x[2])]
  # (div (grad u))(x) = -4*(pi^2)*sin(sqrt(2)*pi*x[1])*cos(sqrt(2)*pi*x[2])
  #f(x::Vector{R}) = 4*(pi^2)*sin(sqrt(2)*pi*x[1])*cos(sqrt(2)*pi*x[2])
  #g(x::Vector{R}) = u(x)
  ############## example 4
  #u(x::Vector{R}) = x[2]*sin(sqrt(2)*pi*x[1])
  #grad_u = [x -> sqrt(2)*pi*x[2]*cos(sqrt(2)*pi*x[1]), x -> sin(sqrt(2)*pi*x[1])]
  # (div (grad u))(x) = -2.*pi*x[2]*sin(sqrt(2)*pi*x[1])
  #f(x::Vector{R}) = 2.*(pi^2)*x[2]*sin(sqrt(2)*pi*x[1])
  #g(x::Vector{R}) = u(x)
############## example 5
  u(x::Vector{R}) = .5*(x[1])^2+sin(sqrt(2)*pi*x[2])
  grad_u = [x -> x[1], x -> sqrt(2)*pi*cos(sqrt(2)*pi*x[2])]
  #(div (grad u))(x) = -2*(pi^2)*sin(sqrt(2)*pi*x[2])+1
  f(x::Vector{R}) = 2*(pi^2)*sin(sqrt(2)*pi*x[2])-1
  g(x::Vector{R}) = u(x)


  const s = 2//1
  const mod_prob = LaplaceModelProblem(f, g, [0.,0.], [1./sqrt(2), 1./sqrt(2)], s, u, grad_u)
  const k = deg(1)
  const r = deg(2)

  const alpha = (k + s - 1)/(r + 1 - min(0, 2 - s))
  const initial_mesh_pair_spec = MeshPairSpec(4,1) # depends on alpha!
  #const initial_mesh_pair_spec = MeshPairSpec(2,1) # depends on alpha!
  const num_mesh_pairs = 3

  # Shouldn't need to change anything from here down.

  println("Expected L2 convergence rate is $((k + s - 1)*(r + 1)/(r + 1 + max(0, s-2))).")
  println("Expected H1 convergence rate is $((k + s - 1)*r/(r + 1 + max(0, s-2))).")

  const mesh_pair_specs = makeMeshPairSpecs(initial_mesh_pair_spec, num_mesh_pairs, alpha)

  test_superconvergence(mod_prob,
                        k,
                        r,
                        mesh_pair_specs)
end
