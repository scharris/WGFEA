module WG
export WGSolver, solve

using Common
import Mesh, Mesh.FENum, Mesh.FERelFace
import VBF, VBF.AbstractVariationalBilinearForm
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.bel_num, WGBasis.MonNum
import Poly.Polynomial

# METHOD
# Let {b_i}_i be a basis for V_h^0(Omega), and vbf the bilinear form for
# the variational problem.  Then the WG approximate solution u_h satisfies
#   vbf(u_h, v) = (f,v_0) for all v in V_h^0(Omega)
# which holds iff
#   vbf(u_h, b_i) = (f, (b_i)_0) for all basis elements b_i
#
# With
#   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on each segment of the outside
# boundary of Omega and 0 elsewhere, this becomes
#   vbf(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i,
# ie.,
#   (sys)
#         sum_j{ vbf(b_j, b_i) eta_j } = (f, (b_i)_0) - vbf(Q_b g, b_i) for all i
#
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.
# Note that the matrix for the system m is given by m_{i,j} = vbf(b_j, b_i).

type WGSolver

  vbf::AbstractVariationalBilinearForm

  vbf_bel_vs_bel_transpose::Matrix{R}

  basis::WeakFunsPolyBasis

  function WGSolver(vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
    new(vbf, VBF.bel_vs_bel_transpose(basis, vbf), basis)
  end
end

function solve(f::Function, g::Function, basis::WeakFunsPolyBasis, wg_solver::WGSolver)
  const rhs = sys_rhs(f, g, wg_solver.vbf)
  wg_solver.vbf_bel_vs_bel_transpose \ rhs
end

# This type is used to cache projection polynomials in a matrix.
PolyOrInt = Union(Polynomial, Int)

# Compute the vector of right hand sides of (sys) for all basis indexes i,
# with the i^th component being (f, (b_i)_0) - vbf(Q_b g, b_i).
function sys_rhs(f::Function, g::Function, vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
  const rhs = Array(R, basis.total_bels)
  const g_projs_cache = zeros(PolyOrInt, Mesh.num_fes(basis.mesh), Mesh.num_side_faces_per_fe(basis.mesh))
  for i=1:basis.total_bels
    const bel_i = bel_num(i)
    rhs[i] = ip_on_interiors(f, bel_i, basis) - vbf_proj_g_on_bsides_vs_bel(vbf, g, bel_i, basis, g_projs_cache)
  end
  rhs
end

function ip_on_interiors(f::Function, bel::BElNum, basis::WeakFunsPolyBasis)
  if WGBasis.is_interior_supported(bel, basis)
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_mon = WGBasis.interior_monomial(bel, basis)
    Mesh.integral_global_x_face_rel_on_fe_face(f, bel_mon, bel_fe, Mesh.interior_face, basis.mesh)
  else
    zeroR
  end
end

# Evaluate the variational bilinear form for the projection of the boundary value function g onto
# outside boundary segments vs. the given basis element.  This is the vbf(Q_b g, b_i) term of
# the right hand side of (sys).
function vbf_proj_g_on_bsides_vs_bel(vbf::AbstractVariationalBilinearForm,
                                     g::Function,
                                     bel::BElNum,
                                     basis::WeakFunsPolyBasis,
                                     g_projs_cache::Matrix{PolyOrInt})
  # By linearity of vbf in its first parameter, we can sum the contributions
  # of the projections onto the individual outside boundary sides.  By the
  # locality assumption of our variational bilinear forms (see the VBF module
  # notes), we only need to consider those outside boundary sides which share
  # a finite element with the support of the given basis element.
  bside_contrs = zeroR
  if WGBasis.is_interior_supported(bel)
    # Only any outside boundary sides which are included in the bel's support fe can contribute.
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_monn = WGBasis.interior_monomial_num(bel, basis)
    for s=fe_face(1):fe_face(Mesh.num_side_faces_per_fe(mesh))
      if Mesh.is_boundary_side(fe, s, mesh)
        const proj_g = proj_on_fe_side(g, bel_fe, s, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_int_mon(vbf, bel_fe, proj_g, s, bel_monn, basis)
      end
    end
  else # side supported bel
    # Only any outside boundary sides which are included in one of the including fe's of the bel support can contribute.
    const supp_incls = WGBasis.fe_inclusions_of_side_support(bel, basis)
    const bel_monn = WGBasis.side_monomial_num(bel, basis)
    # Sum contributions from outside boundary sides of the first including fe.
    for s=fe_face(1):fe_face(Mesh.num_side_faces_per_fe(mesh))
      if Mesh.is_boundary_side(supp_incls.fe1, s, mesh)
        const proj_g = proj_on_fe_side(g, supp_incls.fe1, s, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_side_mon(vbf, supp_incls.fe1, proj_g, s, bel_monn, supp_incls.face_in_fe1, basis)
      end
    end
    # Sum contributions from outside boundary sides of the second including fe.
    for s=fe_face(1):fe_face(Mesh.num_side_faces_per_fe(mesh))
      if Mesh.is_boundary_side(supp_incls.fe2, s, mesh)
        const proj_g = proj_on_fe_side(g, supp_incls.fe2, s, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_side_mon(vbf, supp_incls.fe2, proj_g, s, bel_monn, supp_incls.face_in_fe2, basis)
      end
    end
  end
  bside_contrs
end

function vbf_poly_on_fe_bside_vs_int_mon(vbf::AbstractVariationalBilinearForm,
                                         fe::FENum,
                                         p::Polynomial, p_bside_face::FERelFace,
                                         int_monn::MonNum,
                                         basis::WeakFunsPolyBasis)
  p_mon_contrs = zeroR
  for i=1:length(p.mons)
    const p_monn = WGBasis.mon_num_for_mon_on_side_face(p.mons[i], p_bside_face, basis)
    p_mon_contrs += p.coefs[i] * VBF.side_mon_vs_int_mon(fe, p_monn, p_bside_face, int_monn, basis, vbf)
  end
  p_mon_contrs
end

function vbf_poly_on_fe_bside_vs_side_mon(vbf::AbstractVariationalBilinearForm,
                                          fe::FENum,
                                          p::Polynomial, p_bside_face::FERelFace,
                                          side_monn::MonNum, side_mon_face::FERelFace,
                                          basis::WeakFunsPolyBasis)
  p_mon_contrs = zeroR
  for i=1:length(p.mons)
    const p_monn = WGBasis.mon_num_for_mon_on_side_face(p.mons[i], p_bside_face, basis)
    p_mon_contrs += p.coefs[i] * VBF.side_mon_vs_side_mon(fe, p_monn, p_bside_face, side_monn, side_mon_face, basis, vbf)
  end
  p_mon_contrs
end


# Retrieve the indicated projection from cache or compute it, writing to the cache if it was not already cached.
function proj_on_fe_side(g::Function, fe::FENum, side_face::FERelFace, basis::WeakFunsPolyBasis, g_projs_cache::Matrix{PolyOrInt})
  const cached_proj = g_projs_cache[fe, side_face]
  if cached_proj != 0
    cached_proj
  else
    const proj = Proj.project_onto_fe_face(g, fe, side_face, basis)
    g_projs_cache[fe, side_face] = proj
    proj
  end
end


end # end of module
