module WG
export WGSolver, solve

using Common
import Mesh, Mesh.FENum
import VBF, VBF.AbstractVariationalBilinearForm
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.bel_num
import Poly.PolynomialVector

# METHOD
# Let {b_i}_i be a basis for V_h^0(Omega), and bf the bilinear form for
# the variational problem.  Then the WG approximate solution u_h satisfies
#   bf(u_h, v) = (f,v_0) for all v in V_h^0(Omega)
# which holds iff
#   bf(u_h, b_i) = (f, (b_i)_0) for all basis elements b_i
#
# With
#   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on each segment of the outside
# boundary of Omega and 0 elsewhere, this becomes
#   bf(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i, ie.
#
#   (sys)
#         sum_j{bf(b_j, b_i) eta_j} = (f, (b_i)_0) - bf(Q_b g, b_i) for all i
#
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.
# Note that the matrix for the system m is given by m_{i,j} = bf(b_j, b_i).

type WGSolver

  bilinear_form::AbstractVariationalBilinearForm

  basis::WeakFunsPolyBasis

  bf_transpose_bel_vs_bel::Matrix{R}

  function WGSolver(bf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
    new(bf, basis, VBF.transpose_bel_vs_bel_matrix(basis, bf))
  end
end

function solve(f::Function, g::Function, basis::WeakFunsPolyBasis, wg_solver::WGSolver)
  const rhs = sys_rhs(f, g, wg_solver.bilinear_form)
  wg_solver.bf_transpose_bel_vs_bel \ rhs
end

# This type is used to cache polynomial vectors in a matrix.
PolyVecOrInt = Union(PolynomialVector, Int)

# Compute the vector of right hand sides of (sys) for all basis indexes i, with i^th component being
# (f, (b_i)_0) - bf(Q_b g, b_i).
function sys_rhs(f::Function, g::Function, bf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
  const rhs = Array(R, basis.total_bels)
  const projs_cache = zeros(PolyVecOrInt, Mesh.num_fes(basis.mesh), Mesh.num_side_faces_per_fe(basis.mesh))
  for i=1:basis.total_bels
    const bel = bel_num(i)
    rhs[i] = ip_on_interiors(f, bel, basis) - bf_proj_g_on_b_sides_vs_bel(bf, g, bel, basis, projs_cach)
  end
  rhs
end

function ip_on_interiors(f::Function, bel::BElNum, basis::WeakFunsPolyBasis)
  if WGBasis.is_interior_supported(bel, basis)
    const mon = WGBasis.interior_monomial(bel, basis)
    const fe = WGBasis.support_interior_num(bel, basis)
    Mesh.integral_prod_on_fe_face(f, mon, fe, Mesh.interior_face, basis.mesh)
  else
    zeroR
  end
end

function bf_proj_g_on_b_sides_vs_bel(bf::AbstractVariationalBilinearForm,
                                     g::Function,
                                     bel::BElNum,
                                     basis::WeakFunsPolyBasis,
                                     projs_cache::Matrix{PolyVecOrInt})
  if WGBasis.is_interior_supported(bel)
    const fe = WGBasis.support_interior_num(bel, basis)
    bf_proj_g_on_fe_b_sides_vs_bel(bf, g, fe, bel, basis, projs_cache)
  else # side supported bel
    const nb_side_num = WGBasis.support_nb_side_num(bel, basis)
    const side_incls = fe_inclusions_of_nb_side(nb_side_num, basis.mesh)
    bf_proj_g_on_fe_b_sides_vs_bel(bf, g, side_incls.fe1, bel, basis, projs_cache) +
    bf_proj_g_on_fe_b_sides_vs_bel(bf, g, side_incls.fe2, bel, basis, projs_cache)
  end
end

function bf_proj_g_on_fe_b_sides_vs_bel(bf::AbstractVariationalBilinearForm,
                                        g::Function,
                                        fe::FENum,
                                        bel::BElNum,
                                        basis::WeakFunsPolyBasis,
                                        projs_cache::Matrix{PolyVecOrInt})
  const mesh = basis.mesh
  # Sum the contributions from the individual boundary side projections, which we can do
  # because of the linearity of the bf in its first parameter.
  bside_contrs_sum = zeroR
  for j=1:Mesh.num_side_faces_per_fe(mesh)
    const side = fe_face(j)
    if Mesh.is_boundary_side(fe, side, mesh)
      const proj_g = begin
        const cached_proj = projs_cache[fe, side]
        if cached_proj != 0
          cached_proj
        else
          const proj = Proj.project_onto_fe_face(g, fe, side, basis)
          cached_projs[fe, side] = proj
          proj
        end
      end
      bside_contrs_sum += VBF.poly_on_fe_face_vs_bel(proj_g, fe, side, bel, basis, bf)
    end
  end
  bside_contrs_sum
end

end # end of module
