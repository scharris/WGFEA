module Proj
export project_onto_fe_face, project_local_poly_onto_face

using Common
import Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FEFace
import WGBasis, WGBasis.WeakFunsPolyBasis

# Projection onto a Subspace
# --------------------------
# <Proj(f),e_i> = <f,e_i> for all e_i in the passed subspace basis.
# Since Proj(f) = sum_j{ a_j e_j } for some a_j's, this gives us a linear system
#   sum_j { <e_j,e_i> a_j } = <f, e_i>
# which we can solve for the a_j.

function project_onto_fe_face(g::Function, fe::FENum, face::FEFace, basis::WeakFunsPolyBasis)
  if face == Mesh.interior_face
    const ips_bel_vs_bel = WGBasis.ref_interior_mons_L2ips_matrix(basis)
    const ips_g_vs_bels = [Mesh.integral_prod_on_fe_face(g, mon, fe, face, basis.mesh) for mon in basis.ref_interior_mons]
    const coefs = ips_bel_vs_bel \ ips_g_vs_bels
    Polynomial(basis.ref_interior_mons, coefs)
  else
    const ips_bel_vs_bel = WGBasis.ref_side_mons_L2ips_matrix(basis)
    const ips_g_vs_bels = [Mesh.integral_prod_on_fe_face(g, mon, fe, face, basis.mesh) for mon in basis.ref_side_mons]
    const coefs = ips_bel_vs_bel \ ips_g_vs_bels
    Polynomial(basis.ref_side_mons, coefs)
  end
end

function project_local_poly_onto_face(p::Polynomial, face::FEFace, basis::WeakFunsPolyBasis)
  if face == Mesh.interior_face
    const ips_bel_vs_bel = WGBasis.ref_interior_mons_L2ips_matrix(basis)
    const ips_g_vs_bels = [Mesh.integral_prod_on_ref_fe_face(p, mon, face, basis.mesh) for mon in basis.ref_interior_mons]
    const coefs = ips_bel_vs_bel \ ips_g_vs_bels
    Polynomial(basis.ref_interior_mons, coefs)
  else
    const ips_bel_vs_bel = WGBasis.ref_side_mons_L2ips_matrix(basis)
    const ips_g_vs_bels = [Mesh.integral_prod_on_ref_fe_face(p, mon, face, basis.mesh) for mon in basis.ref_side_mons]
    const coefs = ips_bel_vs_bel \ ips_g_vs_bels
    Polynomial(basis.ref_side_mons, coefs)
  end
end

end # end of module
