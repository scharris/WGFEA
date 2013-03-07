module WGMethod

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
#   u_h = sum_j{eta_j b_j} + Q_b g
# this becomes
#   a_s(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0) for all i
# ie.
#   sum_j{a_s(b_j, b_i) eta_j} = (f, (b_i)_0) - a_s(Q_b g, b_i) for all i
# which is a linear system we can solve for the unknown eta_j coefficients defining u_h.

function wg(bilinear_form::Function, boundary_fn::Function, approx_polys_interiors_deg::Pwr, mesh::RectMesh)
  const k = approx_polys_interiors_deg
  const bf = bilinear_form
  const g = boundary_fn

end # wg method


end # end of module
