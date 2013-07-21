
using Common

import Poly
import RMesh
import TMesh
import WGBasis
import VBF_a_s.a_s
import WGSol


immutable LaplaceModelProblem
  laplace_eq_rhs::Function
  boundary_fn::Function
  region_min::Array{R,1}
  region_max::Array{R,1}
  sol_H_regularity::R
  sol::Function
  grad_sol::Function
end


function test_superconvergence(mod_prob::LaplaceModelProblem,
                               wg_polys_max_deg::Deg,
                               proj_space_polys_max_deg::Deg,
                               proj_mesh_side_secs::Array{Int,1}) # array of # of sections per side for projection mesh
const k = wg_polys_max_deg
const r = proj_space_polys_max_deg
const s = mod_prob.sol_H_regularity
const alpha = (k + s - 1.)/(r + 1. - min(0., 2. - s))

for side_secs in proj_mesh_side_secs

  const tau = norm(mod_prob.region_max - mod_prob.region_min) / side_secs
  const h = tau^(1/alpha)

  const tau_mesh = RectMesh(mod_prob.region_min, mod_prob.region_max, [mesh_coord(num_tau_divs), mesh_coord(num_tau_divs)])
  
  const h_mesh = triangle_submesh(tau_mesh, h)

  const wg_sol = solve_wg(mod_prob, h_mesh, k)


  # TODO: Project the wg sol u_h onto piecewise polynomials of deg <= r on the tau mesh.
  #  Given an element e_tau of the tau mesh, and a basis {v_i} of P_r(e_tau), we have
  #    (Q_tau u_h, v_i)_{e_tau} = (u_h, v_i)_{e_tau} for i=1..dim(P_r(e_tau)).
  #  Letting e_{h,j} be the h-mesh elements whose union is e_tau (by construction),
  #  this gives us
  #    (Q_tau u_h, v_i)_{e_tau} = sum_j (u_h, v_i)_{e_{h,j}}.
  #  Letting lambda be the vector of unknown coefficients satisfying
  #    Q_tau u_h = sum_j lambda_j v_j,
  #  we obtain the linear system
  #    M lambda = R  (sys)
  #      where M_ij = (v_j, v_i)_{e_tau},
  #        and R_i = sum_j (u_h, v_i)_{e_{h,j}} for i=1..dim(P_r(e_tau)).
  const projs_by_tau_fenum = Array(Array{R,1}, num_tau_fes)
  for l=1:num tau mesh els
    sys_rhs = zeros(R, num tau mesh basis mons)
    for each basis mon (# i) of the tau fe
      for each h-mesh fe in e_tau
        const tau_mon_i = ...
        const tau_mon_i_h_rel = ...
        const u_h = WGSol.(...)
        sys_rhs[i] += integral_fe_rel { tau_mon_i_h_rel * u_h }
      end
    end
    projs_by_tau_fenum[l] = # (v_i vs v_j matrix) \ sys_rhs
  end



function triangle_submesh(rmesh::RectMesh, tri_diam::R)
  const num_rmesh_divs = RMesh.mesh_logical_dimensions(rmesh)[1]
  
  const rect_mesh_file = "rect_mesh_$num_rmesh_divs.geo"
  const tri_mesh_file = "tri_mesh_$num_rmesh_divs.msh"

  open(rect_mesh_file, "w") do os
    RMesh.exportAsGmshSurface(os, rmesh, tri_diam)
  end

  run(`gmsh -saveall -2 -clmax $tri_diam -o $tri_mesh_file $rect_mesh_file`)

  return let is = open(tri_mesh_file)
    try 
      TriMesh(is,
              0, # no subdivision
              1e-12, 1e-12, # integral error bounds
              false, # no capture of phys region tags
              true)  # capture geom entity tags (rect mesh fe #'s)
    finally
      close(is)
    end
  end
end

function solve_wg(mod_prob::LaplaceModelProblem, mesh::AbstractMesh, k::Deg)
  const basis = WeakFunsPolyBasis(k, deg(k-1), mesh)
  printBasisSummary(basis)

  println("Constructing VBF (computing vbf bel vs bel matrix)...")
  const vbf = a_s(basis)
  println("VBF constructed.")
  const wg = WGSolver(vbf, basis)

  println("Solving...")
  const wg_sol = solve(mod_prob.laplace_eq_rhs, mod_prob.boundary_fn, wg)
  println("System solved.")

  wg_sol
end

function printBasisSummary(basis::WeakFunsPolyBasis)
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  println("Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements.")
  println("Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println("Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")
  println("Data arrays for non-zeros are of initial (upper bound) size $(int64(interacting_bel_pairs_ub)).")
end