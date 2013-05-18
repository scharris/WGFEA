module WG
export WGSolver, solve

using Common
import Poly, Poly.Polynomial
import Mesh, Mesh.AbstractMesh, Mesh.FENum, Mesh.OrientedShape, Mesh.RelFace, Mesh.rface, Mesh.fen
import Proj
import VBF, VBF.AbstractVariationalBilinearForm
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.beln, WGBasis.MonNum, WGBasis.monn
import WGSol.WGSolution

# METHOD
# Let {b_i}_i be a basis for V_h^0(Omega), and vbf the bilinear form for
# the variational problem.  Then the WG approximate solution u_h satisfies
#   vbf(u_h, v) = (f, v_0) for all v in V_h^0(Omega)
# Since the {b_i}_i are a basis for V_H^0, this holds iff
#   vbf(u_h, b_i) = (f, (b_i)_0) for all basis elements b_i
#
# With
#   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on each segment of
# the outside boundary of Omega and 0 elsewhere, this becomes
#   vbf(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i,
# ie.,
#   (sys)
#         sum_j{ vbf(b_j, b_i) eta_j } = (f, (b_i)_0) - vbf(Q_b g, b_i) for all i
#
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.
# Note that the matrix for the system m is given by m_{i,j} = vbf(b_j, b_i).


immutable WGSolver

  vbf::AbstractVariationalBilinearForm

  vbf_bel_vs_bel_transpose::AbstractMatrix

  basis::WeakFunsPolyBasis

  function WGSolver(vbf::AbstractVariationalBilinearForm, basis::WeakFunsPolyBasis)
    new(vbf, VBF.bel_vs_bel_transpose(basis, vbf), basis)
  end
end

typealias BoundaryProjections Dict{(FENum,RelFace), Vector{R}}


# Solve the system, returning coefficients for all basis elements.
function solve(f::Function, g::FunctionOrConst, wg_solver::WGSolver)
  const g_projs = boundary_projections(g, wg_solver.basis)
  const m = wg_solver.vbf_bel_vs_bel_transpose
  const sol_basis_coefs = m \ sys_rhs(f, g_projs, wg_solver.vbf, wg_solver.basis)
  WGSolution(vec(sol_basis_coefs), g_projs, wg_solver.basis)
end


# Compute the vector of right hand sides of (sys) for all basis indexes i,
# with the i^th component being (f, (b_i)_0) - vbf(Q_b g, b_i).
function sys_rhs(f::Function,
                 g_projs::BoundaryProjections,
                 vbf::AbstractVariationalBilinearForm,
                 basis::WeakFunsPolyBasis)
  const rhs = Array(R, basis.total_bels)
  for i=1:basis.total_bels
    const bel_i = beln(i)
    rhs[i] = ip_on_interiors(f, bel_i, basis) - vbf_boundary_projs_vs_bel(vbf,
                                                                          g_projs,
                                                                          bel_i,
                                                                          basis)
  end
  rhs
end


function ip_on_interiors(f::Function,
                         bel::BElNum,
                         basis::WeakFunsPolyBasis)
  if !WGBasis.is_interior_supported(bel, basis)
    zeroR
  else
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_mon = WGBasis.interior_mon(bel, basis)
    Mesh.integral_global_x_face_rel_on_fe_face(f,
                                               bel_mon,
                                               bel_fe, Mesh.interior_face,
                                               basis.mesh)
  end
end

# Evaluate the variational bilinear form for the projection of the boundary
# value function g onto outside boundary segments vs. the given basis element.
# This is the vbf(Q_b g, b_i) term of the right hand side of (sys). The
# implmentation uses the Element Summability and Locality properties of
# supported variational forms (see VBF module for a discussion).
function vbf_boundary_projs_vs_bel(vbf::AbstractVariationalBilinearForm,
                                   b_projs::BoundaryProjections,
                                   bel::BElNum,
                                   basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  bside_contrs = zeroR
  if WGBasis.is_interior_supported(bel, basis)
    # Only outside boundary sides which are included in the bel's support fe can contribute.
    const bel_fe = WGBasis.support_interior_num(bel, basis)
    const bel_fe_oshape = Mesh.oriented_shape_for_fe(bel_fe, basis.mesh)
    const bel_monn = WGBasis.interior_mon_num(bel, basis)
    for sf=rface(1):Mesh.num_side_faces_for_shape(bel_fe_oshape, mesh)
      if Mesh.is_boundary_side(bel_fe, sf, mesh)
        const proj = boundary_proj(bel_fe, sf, b_projs)
        bside_contrs += vbf_proj_on_oshape_bside_vs_int_mon(vbf,
                                                            proj, bel_fe_oshape, sf,
                                                            bel_monn,
                                                            basis)
      end
    end
  else # side supported bel
    # Only outside boundary sides in one of the fe's including the bel side support can contribute.
    const bel_monn = WGBasis.side_mon_num(bel, basis)
    const bel_supp_incls = WGBasis.fe_inclusions_of_side_support(bel, basis)
    for (bel_supp_fe, bel_supp_rface) in ((bel_supp_incls.fe1, bel_supp_incls.face_in_fe1),
                                          (bel_supp_incls.fe2, bel_supp_incls.face_in_fe2))
      # Sum contributions from outside boundary sides of this supporting fe face.
      for bsf=rface(1):Mesh.num_side_faces_for_fe(bel_supp_fe, mesh) if Mesh.is_boundary_side(bel_supp_fe, bsf, mesh)
        const proj = boundary_proj(bel_supp_fe, bsf, b_projs)
        bside_contrs += vbf_proj_on_fe_bside_vs_side_mon(vbf,
                                                         proj, bel_supp_fe, bsf,
                                                         bel_monn, bel_supp_rface,
                                                         basis)
      end end
    end
  end
  bside_contrs
end

function vbf_proj_on_oshape_bside_vs_int_mon(vbf::AbstractVariationalBilinearForm,
                                             proj_coefs::Vector{R}, fe_oshape::OrientedShape, proj_bside_face::RelFace,
                                             int_monn::MonNum,
                                             basis::WeakFunsPolyBasis)
  proj_mon_contrs = zeroR
  for i=1:length(proj_coefs)
    proj_mon_contrs += proj_coefs[i] * VBF.side_mon_vs_int_mon(fe_oshape,
                                                               monn(i), proj_bside_face,
                                                               int_monn,
                                                               basis,
                                                               vbf)
  end
  proj_mon_contrs
end

function vbf_proj_on_fe_bside_vs_side_mon(vbf::AbstractVariationalBilinearForm,
                                          proj_coefs::Vector{R}, fe::FENum, proj_bside_face::RelFace,
                                          side_monn::MonNum, side_mon_face::RelFace,
                                          basis::WeakFunsPolyBasis)
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, basis.mesh)
  proj_mon_contrs = zeroR
  for i=1:length(proj_coefs)
    proj_mon_contrs += proj_coefs[i] * VBF.side_mon_vs_side_mon(fe_oshape,
                                                                monn(i), proj_bside_face,
                                                                side_monn, side_mon_face,
                                                                basis,
                                                                vbf)
  end
  proj_mon_contrs
end

boundary_proj(fe::FENum, side_face::RelFace, projs::BoundaryProjections) =
  projs[(fe,side_face)]

function boundary_projections(g::FunctionOrConst, basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const projs = Dict{(FENum,RelFace), Vector{R}}()
  sizehint(projs, Mesh.num_boundary_sides(mesh))

  for fe=fen(1):Mesh.num_fes(mesh),
      sf=rface(1):Mesh.num_side_faces_for_fe(fe, mesh)
    if Mesh.is_boundary_side(fe, sf, mesh)
      projs[(fe,sf)] = Proj.project_onto_fe_face(g, fe, sf, basis)
    end
  end

  projs
end

end # end of module
