module WG
export WGSolver, solve

using Common
import Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.fe_num, Mesh.FERelFace, Mesh.fe_face
import Proj
import VBF, VBF.AbstractVariationalBilinearForm
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.bel_num, WGBasis.MonNum

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

# Solve the system, returning all interior and side basis element coefficients, interior coefficients first.
function solve(f::Function, g::Function, wg_solver::WGSolver)
  wg_solver.vbf_bel_vs_bel_transpose \ sys_rhs(f, g, wg_solver.vbf, wg_solver.basis)
end

# Solve the system, returning only interior basis element coefficients.
#function solve(f::Function, g::Function, wg_solver::WGSolver)
#  const all_coefs = solve_full(f, g, wg_solver)
#  all_coefs[1:wg_solver.basis.num_interior_bels]
#end

# This type is used to cache projection polynomials in a matrix.
PolyOrInt = Union(Polynomial, Int)

# Compute the vector of right hand sides of (sys) for all basis indexes i,
# with the i^th component being (f, (b_i)_0) - vbf(Q_b g, b_i).
function sys_rhs(f::Function, g::Function, vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
  const rhs = Array(R, basis.total_bels)
  const g_projs_cache = make_projs_cache(basis)
  for i=1:basis.total_bels
    const bel_i = bel_num(i)
    rhs[i] = ip_on_interiors(f, bel_i, basis) - vbf_proj_g_on_bsides_vs_bel(vbf, g, bel_i, basis, g_projs_cache)
  end
  rhs
end

function ip_on_interiors(f::Function, bel::BElNum, basis::WeakFunsPolyBasis)
  if !WGBasis.is_interior_supported(bel, basis)
    zeroR
  else
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_mon = WGBasis.interior_mon(bel, basis)
    Mesh.integral_global_x_face_rel_on_fe_face(f, bel_mon, bel_fe, Mesh.interior_face, basis.mesh)
  end
end

# Evaluate the variational bilinear form for the projection of the boundary
# value function g onto outside boundary segments vs. the given basis element.
# This is the vbf(Q_b g, b_i) term of the right hand side of (sys). The
# implmentation uses the Element Summability and Locality properties of
# supported variational forms (see VBF module for a discussion).
function vbf_proj_g_on_bsides_vs_bel(vbf::AbstractVariationalBilinearForm,
                                     g::Function,
                                     bel::BElNum,
                                     basis::WeakFunsPolyBasis,
                                     g_projs_cache::Array{Array{PolyOrInt,1}})
  const mesh = basis.mesh
  bside_contrs = zeroR
  if WGBasis.is_interior_supported(bel, basis)
    # Only any outside boundary sides which are included in the bel's support fe can contribute.
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_monn = WGBasis.interior_mon_num(bel, basis)
    for sf=fe_face(1):Mesh.num_side_faces_for_fe(bel_fe, mesh)
      if Mesh.is_boundary_side(bel_fe, sf, mesh)
        const proj_g = proj_on_fe_side(g, bel_fe, sf, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_int_mon(vbf, proj_g, bel_fe, sf, bel_monn, basis)
      end
    end
  else # side supported bel
    # Only outside boundary sides which are included in one of the including fe's of the bel side support can contribute.
    const supp_incls = WGBasis.fe_inclusions_of_side_support(bel, basis)
    const bel_monn = WGBasis.side_mon_num(bel, basis)
    # Sum contributions from outside boundary sides of the first including fe.
    for sf=fe_face(1):Mesh.num_side_faces_for_fe(supp_incls.fe1, mesh)
      if Mesh.is_boundary_side(supp_incls.fe1, sf, mesh)
        const proj_g = proj_on_fe_side(g, supp_incls.fe1, sf, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_side_mon(vbf, proj_g, supp_incls.fe1, sf, bel_monn, supp_incls.face_in_fe1, basis)
      end
    end
    # Sum contributions from outside boundary sides of the second including fe.
    for sf=fe_face(1):Mesh.num_side_faces_for_fe(supp_incls.fe2, mesh)
      if Mesh.is_boundary_side(supp_incls.fe2, sf, mesh)
        const proj_g = proj_on_fe_side(g, supp_incls.fe2, sf, basis, g_projs_cache)
        bside_contrs += vbf_poly_on_fe_bside_vs_side_mon(vbf, proj_g, supp_incls.fe2, sf, bel_monn, supp_incls.face_in_fe2, basis)
      end
    end
  end
  bside_contrs
end

function vbf_poly_on_fe_bside_vs_int_mon(vbf::AbstractVariationalBilinearForm,
                                         p::Polynomial,
                                         fe::FENum,
                                         bside_face::FERelFace,
                                         int_monn::MonNum,
                                         basis::WeakFunsPolyBasis)
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
  p_mon_contrs = zeroR
  const cache = Dict{Any,Any}()
  for i=1:length(p.mons)
    const p_monn = WGBasis.mon_num_for_mon_on_oshape_side(p.mons[i], fe_oshape, bside_face, basis)
    p_mon_contrs += p.coefs[i] * VBF.side_mon_vs_int_mon(fe, p_monn, bside_face, int_monn, basis, cache, vbf)
  end
  p_mon_contrs
end

function vbf_poly_on_fe_bside_vs_side_mon(vbf::AbstractVariationalBilinearForm,
                                          p::Polynomial,
                                          fe::FENum,
                                          p_bside_face::FERelFace,
                                          side_monn::MonNum, side_mon_face::FERelFace,
                                          basis::WeakFunsPolyBasis)
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
  p_mon_contrs = zeroR
  const cache = Dict{Any,Any}()
  for i=1:length(p.mons)
    const p_monn = WGBasis.mon_num_for_mon_on_oshape_side(p.mons[i], fe_oshape, p_bside_face, basis)
    p_mon_contrs += p.coefs[i] * VBF.side_mon_vs_side_mon(fe, p_monn, p_bside_face, side_monn, side_mon_face, basis, cache, vbf)
  end
  p_mon_contrs
end


# Retrieve the indicated projection from cache or compute it, writing to the cache if it was not already cached.
function proj_on_fe_side(g::Function, fe::FENum, side_face::FERelFace, basis::WeakFunsPolyBasis, g_projs_cache::Array{Array{PolyOrInt,1}})
  const cached_proj = g_projs_cache[fe][side_face]
  if cached_proj != 0
    cached_proj
  else
    const proj = Proj.project_onto_fe_face(g, fe, side_face, basis)
    g_projs_cache[fe][side_face] = proj
    proj
  end
end

function make_projs_cache(basis::WeakFunsPolyBasis)
  const num_fes = Mesh.num_fes(basis.mesh)
  const projs_by_fe = Array(Array{PolyOrInt,1}, num_fes)
  for fe=fe_num(1):num_fes
    const num_sides = Mesh.num_side_faces_for_fe(fe, basis.mesh)
    projs_by_fe[fe] = zeros(PolyOrInt, num_sides)
  end
  projs_by_fe
end

end # end of module
