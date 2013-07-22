
using Common

import Poly
import Mesh, Mesh.FENum, Mesh.fenum
import RMesh
import TMesh
import WGBasis
import VBF_a_s.a_s
import WG, WG.WGSolver
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

    const wg_sol = wg_solution(mod_prob, h_mesh, k)

    const projs_by_tau_fe = project(wg_sol, h_mesh, tau_mesh, r)

  end

  # TODO: evaluate the error of the exact solution vs. the projections
end

#  Project the wg sol onto the space of piecewise polynomials of deg <= r on the tau mesh.
#  Returns an array of projection coefficient arrays, indexed by tau mesh element number.
function project(wg_sol::WGSol, h_mesh::AbstractMesh, tau_mesh::RectMesh, r::Deg)
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
  const num_tau_mesh_fes = Mesh.num_fes(tau_mesh)
  const projs_by_tau_fenum = sizehint(Array(Array{R,1},0), num_tau_mesh_fes) # return value

  # Reference basis monomials, which when interpreted locally on a tau mesh element will form
  # a basis of the space to which we are projecting for the element.
  const proj_space_basis_mons = Poly.mons_of_deg_le(r, Mesh.space_dim(tau_mesh))                   # {v_*}
  const proj_space_dim = length(proj_space_basis_mons)

  const included_h_fes_by_tau_fe = fenums_by_geoment(h_mesh, num_tau_mesh_fes)

  const proj_space_basis_vs_basis_ips = ips_matrix(proj_space_basis_mons, RMesh.fe_dims(tau_mesh)) # M
  const sys_rhs = Array(R, proj_space_dim)                                                         # b

  for tau_fe=fenum(1):num_tau_mesh_fes
    fill!(sys_rhs, zeroR)
    # Compute rhs (b) column vector of (sys).
    for i=1:proj_space_dim # (sys) equation number, index into proj space basis
      # Compute component i of (sys) rhs.
      const tau_mon = proj_space_basis_mons[i]                                                     # v_i
      const h_fes = included_h_fes_by_tau_fe[tau_fe]                                               # e_{h,*}
      for h_fe in h_fes                                                                            # e_{h,j}
        const u_h_rel = WGSol.(...) TODO: get poly representing wg sol on h_fe
        # Translate the tau monomial to a polynomial to be interpreted relative to the h_fe origin.
        const tau_mon_rel = Poly.translate(tau_mon, origin_h_fe - origin_tau_fe)
        sys_rhs[i] += Mesh.integral_face_rel_on_oshape_face(u_h_rel * tau_mon_rel, Mesh.oriented_shape_for_fe(h_fe), Mesh.interior_face)
      end
    end
    const proj_coefs = proj_space_basis_vs_basis_ips \ sys_rhs
    push!(projs_by_tau_fenum, proj_coefs)
  end

  assert(length(projs_by_tau_fenum) == num_tau_mesh_fes)
  projs_by_tau_fenum
end

function triangle_submesh(rmesh::RectMesh, tri_diam::R)
  const num_rmesh_divs = RMesh.mesh_ldims(rmesh)[1]
  
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

function wg_solution(mod_prob::LaplaceModelProblem, mesh::AbstractMesh, k::Deg)
  const basis = WGBasis.WeakFunsPolyBasis(k, deg(k-1), mesh)
  printBasisSummary(basis)

  println("Constructing VBF (computing vbf bel vs bel matrix)...")
  const vbf = a_s(basis)
  println("VBF constructed.")
  const wg = WG.WGSolver(vbf, basis)

  println("Solving...")
  const wg_sol = WG.solve(mod_prob.laplace_eq_rhs, mod_prob.boundary_fn, wg)
  println("System solved.")

  wg_sol
end

function fenums_by_geoment(tri_mesh::TMesh, num_geoments::FENum)
  const fenums_by_tag = sizehint(Dict{TMesh.Tag, Array{FENum,1}}(), num_geoments)
  const approx_fes_per_geoment = int64(Mesh.num_fes(tri_mesh) / num_geoments)
  for fe=fenum(1):Mesh.num_fes(tri_mesh)
    const ge = TMesh.geometric_entity_tag(fe, tri_mesh)
    ge_fes = get(d, ge, nothing)
    if ge_fes == nothing
      ge_fes = Array(FENum, int(approx_fes_per_geoment * 1.10))
      d[ge] = ge_fes
    end
    push!(ge_fes, fe)
  end
  d
end

function ips_matrix(mons::Array{Monomial,1}, fe_dims::Array{R,1})
  #TODO
end

function printBasisSummary(basis::WeakFunsPolyBasis)
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  println("Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements.")
  println("Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println("Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")
  println("Data arrays for non-zeros are of initial (upper bound) size $(int64(interacting_bel_pairs_ub)).")
end