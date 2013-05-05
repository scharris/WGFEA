module WGSol
export WGSolution,
       wg_sol_at_interior_rel_point,
       wg_sol_interior_poly,
       wg_sol_side_poly,
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
import VBF, VBF.AbstractVariationalBilinearForm

immutable WGSolution
  basis_coefs::Vector{R}
  boundary_projs::Dict{(FENum,FERelFace), Vector{R}}
  basis::WeakFunsPolyBasis
end


# Solution evaluation functions

function wg_sol_at_interior_rel_point(fe::FENum, x::Vector{R}, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const int_mons = WGBasis.interior_mons(basis)
  sum = zeroR
  for monn=mon_num(1):WGBasis.mons_per_fe_interior(basis)
    const beln = WGBasis.interior_mon_bel_num(fe, monn, basis)
    sum += sol_coefs[beln] * Poly.value_at(int_mons[monn], x)
  end
  sum
end

# Return interior-relative polynomial coefficients which represent the passed WG solution on the
# indicated finite element interior.
wg_sol_interior_coefs(fe::FENum, wg_sol::WGSolution) =
  WGBasis.fe_interior_poly_coefs(fe, wg_sol.basis_coefs, wg_sol.basis)

wg_sol_interior_poly(fe::FENum, wg_sol::WGSolution) =
  Polynomial(WGBasis.interior_mons(wg_sol.basis),
             WGBasis.fe_interior_poly_coefs(fe, wg_sol.basis_coefs, wg_sol.basis))

wg_sol_side_poly(fe::FENum, side_face::FERelFace, wg_sol::WGSolution) =
  Polynomial(WGBasis.side_mons_for_fe_side(fe, side_face, wg_sol.basis),
             WGBasis.fe_side_poly_coefs(fe, side_face, wg_sol.basis_coefs, wg_sol.basis))

wg_sol_side_coefs(fe::FENum, side_face::FERelFace, wg_sol::WGSolution) =
  WGBasis.fe_side_poly_coefs(fe, side_face, wg_sol.basis_coefs, wg_sol.basis)


# Return an interior-relative polynomial vector which represents the weak gradient of the passed
# solution on the indicated finite element interior.
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

function err_L2_norm(exact_sol::Function, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const one = Mesh.one_mon(mesh)
  const d = Mesh.space_dim(mesh)
  const fe_rel_x = Array(R, d)
  const fe_origin = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    Mesh.fe_interior_origin!(fe, fe_origin, mesh)
    function sq_diff(x::Vector{R})
      for i=1:d  fe_rel_x[i] = x[i] - fe_origin[i] end
      const diff = exact_sol(x) - wg_sol_at_interior_rel_point(fe, fe_rel_x, wg_sol)
      diff * diff
    end
    sum_fe_diff_norm_sqs += Mesh.integral_global_x_face_rel_on_fe_face(sq_diff, one, fe, Mesh.interior_face, mesh)
  end
  sqrt(sum_fe_diff_norm_sqs)
end


function err_vs_proj_L2_norm(exact_sol::Function, wg_sol::WGSolution)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const fe_rel_x = Array(R, Mesh.space_dim(mesh))

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const proj_poly = Proj.project_onto_fe_face_as_poly(exact_sol, fe, Mesh.interior_face, basis)
    const wg_sol_poly = wg_sol_interior_poly(fe, wg_sol)
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
  const fe_origin = Array(R, d)

  sum_fe_diff_norm_sqs = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    Mesh.fe_interior_origin!(fe, fe_origin, mesh)
    const wgrad = wg_sol_wgrad_on_interior(fe, wg_sol)
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

# TODO: Need to distinguish boundary sides and non-boundary sides.  We can assume that the boundary sides don't
# contribute here, because f = g on the outside boundary and u_h has g projection values there.
# Compute semi-norm ||| Q_h u - u_h |||_vbf = sqrt(vbf(Q u - u_h, Q u - u_h)),
# where Q_h is piecwise L2 projection onto finite element faces.
function err_vs_proj_vbf_seminorm(exact_sol::Function, wg_sol::WGSolution, vbf::AbstractVariationalBilinearForm)
  const basis = wg_sol.basis
  const mesh = basis.mesh
  const vbf_symm = VBF.is_symmetric(vbf)

  # work arrays
  const side_diffs = Array(Vector{R}, Mesh.max_num_shape_sides(basis.mesh))
  const is_nb_side = Array(Bool, Mesh.max_num_shape_sides(basis.mesh))

  # We must evaluate
  #   sum_{fe} vbf(diff(fe), diff(fe)),
  # where
  #   diff(fe) = diff_int(fe) + sum_{s in nb-sides(fe)} diff_s,
  # and diff_f for a face f is the L2 projection onto f of the exact solution minus the wg
  # solution on the face, and 0 elsewhere.
  # By the bilinearity of the vbf, we need to sum the contributions from each finite element
  # of the vbf applied to the differences on each interior or non-boundary side against that
  # on each interior or non-boundary side.

  sum = zeroR
  for fe=fe_num(1):Mesh.num_fes(mesh)
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
    const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, mesh)
    # Find non-boundary sides
    for sf=fe_face(1):num_sides
      is_nb_side[sf] = !Mesh.is_boundary_side(fe, sf, mesh)
    end

    # Compute exact solution projection minus wg solution for all faces.
    const int_diff = Proj.project_onto_fe_face(exact_sol, fe, Mesh.interior_face, basis) - wg_sol_interior_coefs(fe, wg_sol)
    for sf=fe_face(1):num_sides if is_nb_side[sf]
      diff = Proj.project_onto_fe_face(exact_sol, fe, sf, basis)
      diff -= wg_sol_side_coefs(fe, sf, wg_sol)
      side_diffs[sf] = diff
    end end

    # interior vs interior
    sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape, int_diff, Mesh.interior_face, int_diff, Mesh.interior_face, basis, vbf)

    for sf=fe_face(1):num_sides if is_nb_side[sf]
      const sf_diff = side_diffs[sf]

      # Add contributions from pairings of a side and interior.

      # interior vs side
      const int_vs_side = VBF.poly_on_face_vs_poly_on_face(fe_oshape, int_diff, Mesh.interior_face, sf_diff, sf, basis, vbf)
      sum += int_vs_side

      # side vs interior
      sum += if vbf_symm  int_vs_side
             else VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, int_diff, Mesh.interior_face, basis, vbf) end

      # Add contributions for pairings of sides.

      # side vs side
      if vbf_symm
        # contributions from pairings of unequal sides (symmetric case)
        for sf2=fe_face(sf+1):num_sides if is_nb_side[sf2]
          sum += 2 * VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, side_diffs[sf2], sf2, basis, vbf)
        end end
        # contributions from side self-pairings (symmetric case)
        sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, sf_diff, sf, basis, vbf)
      else # vbf not symmetric
        # contributions for all side vs side pairings (asymmetric case)
        for sf2=fe_face(1):num_sides if is_nb_side[sf2]
          sum += VBF.poly_on_face_vs_poly_on_face(fe_oshape, sf_diff, sf, side_diffs[sf2], sf2, basis, vbf)
        end end
      end
    end end # non-boundary sides
  end
  sqrt(sum)
end

end
