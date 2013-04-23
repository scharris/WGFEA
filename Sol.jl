module Sol
export Solution,
       solution_value_at_interior_rel_point,
       solution_poly_on_interior,
       solution_wgrad_on_interior

using Common
import Poly, Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face, Mesh.fe_num
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.mon_num
import VBF.AbstractVariationalBilinearForm

type Solution
  basis_coefs::Vector{R}
  boundary_projs::Dict{(FENum,FERelFace), Vector{R}}
end


# Solution evaluation functions

# For the passed array of solution basis coefficients, which should include coefficients for at least
# all interior basis elements, evaluate the solution at the indicated finite element relative point.
function solution_value_at_interior_rel_point(sol_coefs::Vector{R}, fe::FENum, x::Vector{R}, basis::WeakFunsPolyBasis)
  const int_mons = basis.interior_mons
  const num_int_mons = basis.mons_per_fe_interior
  const first_beln = WGBasis.first_bel_supported_on_fe_interior(fe, basis)
  sum = zeroR
  for i=1:num_int_mons
    sum += sol_coefs[first_beln + i-1] * Poly.value_at(int_mons[i], x)
  end
  sum
end

solution_value_at_interior_rel_point(sol::Solution, fe::FENum, x::Vector{R}, basis::WeakFunsPolyBasis) =
  solution_value_at_interior_rel_point(sol.basis_coefs, fe, x, basis)

# Return an interior-relative polynomial which represents on the indicated finite element interior
# the solution having the passed sequence of basis element coefficients. The basis coefficients
# should include those for at least all interior basis elements.
function solution_poly_on_interior(sol_coefs::Vector{R}, fe::FENum, basis::WeakFunsPolyBasis)
  const int_mons = basis.interior_mons
  const num_int_mons = basis.mons_per_fe_interior
  const first_beln = WGBasis.first_bel_supported_on_fe_interior(fe, basis)
  const coefs = sol_coefs[first_beln : first_beln + num_int_mons - 1]
  Polynomial(int_mons, coefs)
end

solution_poly_on_interior(sol::Solution, fe::FENum, basis::WeakFunsPolyBasis) =
  solution_poly_on_interior(sol.basis_coefs, fe, basis)

# Return an interior-relative polynomial which represents on the indicated finite element interior
# the weak gradient of the solution having the passed sequence of basis element coefficients. The
# basis coefficients should include those for at least all interior basis elements.
function solution_wgrad_on_interior(sol::Solution, fe::FENum, basis::WeakFunsPolyBasis)
  const wgrad_solver = basis.wgrad_solver
  const mesh = basis.mesh
  const sol_coefs = sol.basis_coefs
  const boundary_projs = sol.boundary_projs
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # The weak gradient of the solution is the sum of the weak gradients of the functions which equal
  # the solution on one face of our finite element and are zero elsewhere.
  total_wgrad = Poly.zero_poly_vec(Mesh.space_dim(mesh))

  # Add contributions for the terms of the solution's interior polynomial.
  const first_bel_on_fe_int = WGBasis.first_bel_supported_on_fe_interior(fe, basis)
  for monn=mon_num(1):basis.mons_per_fe_interior
    const beln = first_bel_on_fe_int + monn - 1
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
      const first_bel_on_fe_side = WGBasis.first_bel_supported_on_fe_side(fe, sf, basis)
      for monn=mon_num(1):basis.mons_per_fe_side
        const beln = first_bel_on_fe_side + monn - 1
        total_wgrad += sol_coefs[beln] * WGBasis.wgrad_side_mon(monn, fe_oshape, sf, basis)
      end
    end
  end

  total_wgrad
end


# Error Estimates


function err_est_L2(approx_sol_coefs::Vector{R}, actual_sol::Function, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)

  sum_sq_diff_norms = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const approx_poly = solution_poly_on_interior(approx_sol_coefs, fe, basis)
    const fe_origin = Mesh.fe_interior_origin(fe, mesh)
    function sq_diff(x::Vector{R})
      for i=1:d
        fe_rel_x[i] = x[i] - fe_origin[i]
      end
      const diff = actual_sol(x) - Poly.value_at(approx_poly, fe_rel_x)
      diff * diff
    end
    sum_sq_diff_norms += Mesh.integral_global_x_face_rel_on_fe_face(sq_diff, one, fe, Mesh.interior_face, mesh)
  end
  sqrt(sum_sq_diff_norms)
end

err_est_L2(approx_sol::Solution, actual_sol::Function, basis::WeakFunsPolyBasis) =
  err_est_L2(approx_sol.basis_coefs, actual_sol, basis)

function err_est_vbf(approx_sol::Vector{R}, actual_sol::Function, vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
# TODO
end

end
