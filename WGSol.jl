module WGSol
export WGSolution,
       wg_sol_at_interior_rel_point,
       wg_sol_poly_on_interior,
       wg_sol_wgrad_on_interior,
       err_L2_norm,
       err_wgrad_vs_grad_L2_norm,
       err_vbf_seminorm

using Common
import Poly, Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face, Mesh.fe_num
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.mon_num
import VBF.AbstractVariationalBilinearForm

type WGSolution
  basis_coefs::Vector{R}
  boundary_projs::Dict{(FENum,FERelFace), Vector{R}}
end


# Solution evaluation functions

wg_sol_at_interior_rel_point(wg_sol::WGSolution, fe::FENum, x::Vector{R}, basis::WeakFunsPolyBasis) =
  wg_sol_at_interior_rel_point(wg_sol.basis_coefs, fe, x, basis)

function wg_sol_at_interior_rel_point(sol_coefs::Vector{R}, fe::FENum, x::Vector{R}, basis::WeakFunsPolyBasis)
  const int_mons = basis.interior_mons
  sum = zeroR
  for monn=mon_num(1):basis.mons_per_fe_interior
    const beln = WGBasis.interior_mon_bel_num(fe, monn, basis)
    sum += sol_coefs[beln] * Poly.value_at(int_mons[monn], x)
  end
  sum
end

# Return an interior-relative polynomial which represents on the indicated finite element interior
# the solution having the passed sequence of basis element coefficients. The basis coefficients
# should include those for at least all interior basis elements.
function wg_sol_poly_on_interior(sol_coefs::Vector{R}, fe::FENum, basis::WeakFunsPolyBasis)
  const int_mons = basis.interior_mons
  const num_int_mons = basis.mons_per_fe_interior
  const first_beln = WGBasis.interior_mon_bel_num(fe, mon_num(1), basis)
  # TODO: shouldn't assume this layout of basis elements, defer to WGBasis instead
  const coefs = sol_coefs[first_beln : first_beln + num_int_mons - 1]
  Polynomial(int_mons, coefs)
end

wg_sol_poly_on_interior(wg_sol::WGSolution, fe::FENum, basis::WeakFunsPolyBasis) =
  wg_sol_poly_on_interior(wg_sol.basis_coefs, fe, basis)

# Return an interior-relative polynomial which represents on the indicated finite element interior
# the weak gradient of the solution having the passed sequence of basis element coefficients. The
# basis coefficients should include those for at least all interior basis elements.
function wg_sol_wgrad_on_interior(wg_sol::WGSolution, fe::FENum, basis::WeakFunsPolyBasis)
  const wgrad_solver = basis.wgrad_solver
  const mesh = basis.mesh
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

function err_L2_norm(wg_sol_coefs::Vector{R}, exact_sol::Function, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const wg_sol_poly = wg_sol_poly_on_interior(wg_sol_coefs, fe, basis)
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

err_L2_norm(wg_sol::WGSolution, exact_sol::Function, basis::WeakFunsPolyBasis) =
  err_L2_norm(wg_sol.basis_coefs, exact_sol, basis)

# TODO
function err_wgrad_vs_grad_L2_norm(wg_sol::WGSolution, exact_sol_grad::Function, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const wgrad = wg_sol_wgrad_on_interior(wg_sol, fe, basis)
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

function err_vbf_seminorm(wg_sol::Vector{R}, exact_sol::Function, vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
# TODO
end

end
