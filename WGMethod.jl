module WGMethod
export WGSolver, solve

using Common
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.mon_num
import Mesh, Mesh.AbstractMesh, Mesh.fe_face
import WGrad, WGrad.WGradSolver
import Poly, Poly.PolynomialVector

# METHOD
# Let {b_i}_i be a basis for V_h^0(Omega), and a_s the bilinear form for
# the variational problem.  Then the WG approximate solution u_h satisfies
#   a_s(u_h, v) = (f,v_0) for all v in V_h^0(Omega)
# which holds iff
#   a_s(u_h, b_i) = (f, (b_i)_0) for all b_i
#
# With
#   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on the outside boundary of Omega
# and 0 elsewhere, this becomes
#   a_s(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i
# ie.
#   sum_j{a_s(b_j, b_i) eta_j} = (f, (b_i)_0) - a_s(Q_b g, b_i) for all i
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.

type WGSolver
  bilinear_form::Function
  int_polys_max_deg::Deg
  mesh::AbstractMesh
  basis::WeakFunsPolyBasis
  wgrad_solver::WGradSolver
  int_wgrads_by_mon_num::Array{PolynomialVector,1}
  side_wgrads_by_side_mon_num::Matrix{PolynomialVector}

  function WGSolver(bilinear_form::Function,
                    int_polys_max_deg::Deg,
                    mesh::AbstractMesh)
    const d = Mesh.space_dim(mesh)
    const basis = WeakFunsPolyBasis(int_polys_max_deg, d, mesh)
    const wgrad_solver = WGradSolver(deg(int_polys_max_deg-1), d, mesh)

    # Precompute weak gradients by face and monomial
    const int_wgrads = interior_supported_wgrads_by_mon_num(basis, wgrad_solver)
    const side_wgrads = side_supported_wgrads_by_side_mon_num(basis, wgrad_solver)

    new(bilinear_form,
        int_polys_max_deg,
        mesh,
        basis,
        wgrad_solver,
        int_wgrads,
        side_wgrads)
  end
end


function interior_supported_wgrads_by_mon_num(basis::WeakFunsPolyBasis, wgrad_solver::WGradSolver)
  const wgrads = Array(PolynomialVector, WGBasis.num_monomials_per_fe_interior(basis))
  for m=1:length(wgrads)
    mon = WGBasis.interior_monomial_by_num(mon_num(m), basis)
    wgrads[m] = WGrad.wgrad(mon, Mesh.interior_face, wgrad_solver)
  end
  wgrads
end

function side_supported_wgrads_by_side_mon_num(basis::WeakFunsPolyBasis, wgrad_solver::WGradSolver)
  const sides_per_fe = Mesh.num_side_faces_per_fe(basis.mesh)
  const mons_per_side = WGBasis.num_monomials_per_fe_side(basis)
  const wgrads = Array(PolynomialVector, sides_per_fe, mons_per_side)
  for s=1:sides_per_fe, m=1:mons_per_side
    const mon = WGBasis.side_monomial_by_num(mon_num(m), basis)
    wgrads[s,m] = WGrad.wgrad(mon, fe_face(s), wgrad_solver)
  end
  wgrads
end


function solve(variational_rhs_fun::Function, boundary_fun::Function, solver::WGSolver)
  const f = variational_rhs_fun
  const g = boundary_fun
  const bf = solver.bilinear_form

# TODO: When building main basis vs. basis matrix, loop separately over the interiors and sides,
#       for both i and j, so in each case we know already the type for both and can have the appropriate
#       structures allocated.

#a_s(b_j, b_i)
# - interior vs interior
#    Different interior numbers
#    If b_i and b_j do not have the same supporting interior number, then 0,
#    else
#      The common supporting fe T becomes the context of all ips.
#      Find weak gradient of b_i and b_j on T [via wgrad(b_mon, interior_face, wgrad_solver)].
#      Multiply grad_w b_i with matrix a (new PolynomialVector op).  Can leave this out for now.
#      Compute (a grad_w b_i, grad_w b_j)_T [via Mesh.integral_on_ref_fe_interior(dot_prod, mesh)]
#      Compute stabilization term
#        The values on the actual boundary, (b_i/j)_b are 0, leaving:
#        <b_i, b_j>_boundary(T)
#        [Need new mesh function for this (or new code to represent "all sides"), to integrate a scalar function over all sides of ref fe]

end # wg method


end # end of module
