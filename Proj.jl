module Proj
export project_onto_fe_face,
       project_interior_monomial_onto_side_face

using Common
import Poly, Poly.Polynomial, Poly.Monomial
import Mesh, Mesh.FENum, Mesh.FERelFace
import WGBasis, WGBasis.WeakFunsPolyBasis

# Projection onto a Subspace
# --------------------------
# <Proj(g), e_i> = <g, e_i> for all e_i in the passed subspace basis.
# Since Proj(g) = sum_j{ a_j e_j } for some a_j's, this gives us a linear system
#   sum_j { <e_j, e_i> a_j } = <g, e_i>
# which we can solve for the a_j.

function project_onto_fe_face(g::Function, fe::FENum, face::FERelFace, basis::WeakFunsPolyBasis)
  const proj_poly =
    if face == Mesh.interior_face
      const ref_interior_mons = WGBasis.ref_interior_mons(basis)
      const ips_bel_vs_bel = WGBasis.ips_ref_interior_mons(basis)
      const ips_g_vs_bels = [Mesh.integral_global_x_face_rel_on_fe_face(g, mon, fe, face, basis.mesh)::R
                             for mon in ref_interior_mons]
      const coefs = ips_bel_vs_bel \ ips_g_vs_bels
      Polynomial(ref_interior_mons, coefs)
    else
      const ref_side_mons = WGBasis.ref_side_mons(face, basis)
      const ips_bel_vs_bel = WGBasis.ips_ref_side_mons(face, basis)
      const ips_g_vs_bels = [Mesh.integral_global_x_face_rel_on_fe_face(g, mon, fe, face, basis.mesh)::R
                             for mon in ref_side_mons]
      const coefs = ips_bel_vs_bel \ ips_g_vs_bels
      Polynomial(ref_side_mons, coefs)
    end
  Poly.canonical_form(proj_poly)
end

function project_interior_monomial_onto_side_face(int_mon::Monomial, side_face::FERelFace, basis::WeakFunsPolyBasis)
  const ref_side_mons = WGBasis.ref_side_mons(side_face, basis)
  const ips_side_bels = WGBasis.ips_ref_side_mons(side_face, basis)
  const ips_int_mon_vs_bels = [Mesh.integral_fe_rel_x_side_rel_on_side(int_mon, side_bel_mon, side_face, basis.mesh)::R
                               for side_bel_mon in ref_side_mons]
  const coefs = ips_side_bels \ ips_int_mon_vs_bels
  Polynomial(ref_side_mons, coefs)
end

end # end of module
