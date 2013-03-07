module WG

using WGBasis
using RMesh
using Poly

# a_s(u_h, v) = (f,v_0) for all v in V_h^0(Omega)
# or equivalently,
# a_s(u_h, b_i) = (f, (b_i)_0) for all basis elements b_i
#
# u_h = sum_j{eta_j b_j} + Q_b g
#
# a_s(sum_j{eta_j b_j} + Q_b g, b_i) = (f, (b_i)_0)
#
# sum_j{eta_j a_s(b_j, b_i)} + a_s(Q_b g, b_i) = (f, (b_i)_0)
#
# sum_j{a_s(b_j, b_i) eta_j} = (f, (b_i)_0) - a_s(Q_b g, b_i) for all i
#
# which is a linear system we can solve for the unknowns eta_j coefficients defining u_h.

function wg(f::Function, mesh::RectMesh, k::Pwr, dom_dim::Dim)



end # wg method


end # end of module
