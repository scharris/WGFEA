module WGSol
export WGSolution,
       wg_sol_at_interior_rel_point,
       wg_sol_poly_on_interior,
       wg_sol_wgrad_on_interior,
       err_L2_norm,
       err_vs_proj_L2_norm,
       err_grad_vs_wgrad_L2_norm,
       err_vbf_seminorm

using Common
import Poly, Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face, Mesh.fe_num
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.mon_num
import Proj
import VBF.AbstractVariationalBilinearForm

immutable WGSolution
  basis_coefs::Vector{R}
  boundary_projs::Dict{(FENum,FERelFace), Vector{R}}
  basis::WeakFunsPolyBasis
end


# Solution evaluation functions

wg_sol_at_interior_rel_point(fe::FENum, x::Vector{R}, wg_sol::WGSolution) =
  wg_sol_at_interior_rel_point(fe, x, wg_sol.basis_coefs, wg_sol.basis)

function wg_sol_at_interior_rel_point(fe::FENum, x::Vector{R}, sol_coefs::Vector{R}, basis::WeakFunsPolyBasis)
  const int_mons = basis.interior_mons
  sum = zeroR
  for monn=mon_num(1):basis.mons_per_fe_interior
    const beln = WGBasis.interior_mon_bel_num(fe, monn, basis)
    sum += sol_coefs[beln] * Poly.value_at(int_mons[monn], x)
  end
  sum
end

# Return an interior-relative polynomial which represents on the indicated finite element interior
# the solution having the passed sequence of basis element coefficients.
wg_sol_poly_on_interior(fe::FENum, wg_sol::WGSolution) =
  WGBasis.fe_interior_poly_for_weak_fun_coefs(fe, wg_sol.basis_coefs, wg_sol.basis)

wg_sol_poly_on_interior(fe::FENum, wg_sol_coefs::Vector{R}, basis::WeakFunsPolyBasis) =
  WGBasis.fe_interior_poly_for_weak_fun_coefs(fe, wg_sol_coefs, basis)

# Return an interior-relative polynomial vector which represents on the indicated finite element interior
# the weak gradient of the solution having the passed sequence of basis element coefficients.
function wg_sol_wgrad_on_interior(fe::FENum, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const wgrad_solver = basis.wgrad_solver
  const sol_coefs = wg_sol.basis_coefs
  const boundary_projs = wg_sol.boundary_projs
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # The weak gradient of the solution is the sum of the weak gradients of the functions which equal
  # the solution on one face of our finite element and are zero elsewhere.
  total_wgrad = Poly.zero_poly_vec(Mesh.space_dim(mesh))

  # Add contributions for the terms of the solution's interior polynomial.
  for monn=mon_num(1):basis.mons_per_fe_interior
    const beln = WGBasis.interior_mon_bel_num(fe, monn, basis)
    total_wgrad += sol_coefs[beln] * WGBasis.wgrad_interior_mon(monn, fe_oshape, basis)
  end

  # Add contributions for the solution's side polynomials, for both non-boundary and boundary sides.
  for sf=fe_face(1):Mesh.num_side_faces_for_shape(fe_oshape, mesh)
    if Mesh.is_boundary_side(fe, sf, mesh)
      const proj_coefs = boundary_projs[(fe, sf)]
      for monn=mon_num(1):basis.mons_per_fe_side
        total_wgrad += proj_coefs[monn] * WGBasis.wgrad_side_mon(monn, fe_oshape, sf, basis)
      end
    else
      for monn=mon_num(1):basis.mons_per_fe_side
        const beln = WGBasis.side_mon_bel_num(fe, sf, monn, basis)
        total_wgrad += sol_coefs[beln] * WGBasis.wgrad_side_mon(monn, fe_oshape, sf, basis)
      end
    end
  end

  total_wgrad
end


# Error Estimates

err_L2_norm(exact_sol::Function, wg_sol::WGSolution) =
  err_L2_norm(exact_sol, wg_sol.basis_coefs, wg_sol.basis)

function err_L2_norm(exact_sol::Function, wg_sol_coefs::Vector{R}, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const wg_sol_poly = wg_sol_poly_on_interior(fe, wg_sol_coefs, basis)
    const fe_origin = Mesh.fe_interior_origin(fe, mesh)
    function sq_diff(x::Vector{R})
      for i=1:d  fe_rel_x[i] = x[i] - fe_origin[i] end
      const diff = exact_sol(x) - Poly.value_at(wg_sol_poly, fe_rel_x)
      diff * diff
    end
    sum_fe_diff_norm_sqs += Mesh.integral_global_x_face_rel_on_fe_face(sq_diff, one, fe, Mesh.interior_face, mesh)
  end
  sqrt(sum_fe_diff_norm_sqs)
end


err_vs_proj_L2_norm(exact_sol::Function, wg_sol::WGSolution) =
  err_vs_proj_L2_norm(exact_sol, wg_sol.basis_coefs, wg_sol.basis)

function err_vs_proj_L2_norm(exact_sol::Function, wg_sol_coefs::Vector{R}, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const fe_rel_x = Array(R, Mesh.space_dim(mesh))

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const proj_poly = Proj.project_onto_fe_face_as_poly(exact_sol, fe, Mesh.interior_face, basis)
    const wg_sol_poly = wg_sol_poly_on_interior(fe, wg_sol_coefs, basis)
    const diff = proj_poly - wg_sol_poly
    const diff_sq = diff * diff
    sum_fe_diff_norm_sqs += Mesh.integral_face_rel_on_face(diff_sq, fe, Mesh.interior_face, mesh)
  end
  sqrt(sum_fe_diff_norm_sqs)
end


function err_grad_vs_wgrad_L2_norm(exact_sol_grad::Function, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const wgrad = wg_sol_wgrad_on_interior(fe, wg_sol)
    const fe_origin = Mesh.fe_interior_origin(fe, mesh)
    function diff_norm_sq(x::Vector{R})
      for i=1:d
        fe_rel_x[i] = x[i] - fe_origin[i]
      end
      const diff = exact_sol_grad(x) - Poly.value_at(wgrad, fe_rel_x)
      dot(diff, diff)
    end
    sum_fe_diff_norm_sqs += Mesh.integral_global_x_face_rel_on_fe_face(diff_norm_sq, one, fe, Mesh.interior_face, mesh)
  end
  sqrt(sum_fe_diff_norm_sqs)
end

function err_vs_proj_vbf_seminorm(exact_sol::Function, wg_sol::Vector{R}, vbf::AbstractVariationalBilinearForm)
  # TODO: For each fe,
  #  Project exact sol to fe interior and sides (store latter in an array)
  #  Compute interior difference and side differences as polynomials.
  #
  #  Sum c1 c2 vbf(mon1, mon2) for each pair of terms ci moni from each interior and side pairing.
  # sqrt sum over fes
end

end
