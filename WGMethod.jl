module WGMethod

using Common
using WGBasis
using RMesh
using Poly

# METHOD
# Let {b_i}_i be a basis for V_h^0(Omega), and a_s the bilinear form for
# the variational problem.  Then
#   a_s(u_h, v) = (f,v_0) for all v in V_h^0(Omega)
# iff
#   a_s(u_h, b_i) = (f, (b_i)_0) for all b_i
#
# With
#   u_h = sum_j{eta_j b_j} + Q_b g, where Q_b is L2 projection on the outside boundary of Omega
# and 0 elsewhere, this becomes
#   a_s(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i
# ie.
#   sum_j{a_s(b_j, b_i) eta_j} = (f, (b_i)_0) - a_s(Q_b g, b_i) for all i
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.

function wg(bilinear_form::Function, boundary_fn::Function, approx_polys_interiors_deg::Deg, mesh::RectMesh)
  const k = approx_polys_interiors_deg
  const bf = bilinear_form
  const g = boundary_fn
# TODO: When building main basis vs. basis matrix, loop separately over the interiors and sides,
#       for both i and j, so in each case we know already the type for both and can have the appropriate
#       structures allocated.
end # wg method


end # end of module
