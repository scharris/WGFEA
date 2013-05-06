module Proj
export project_onto_fe_face,
       project_onto_fe_face_as_poly,
       project_interior_mon_onto_oshape_side

using Common
import Poly.Monomial, Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.RelFace, Mesh.OrientedShape
import WGBasis, WGBasis.WeakFunsPolyBasis

# Projection onto a Subspace
# --------------------------
# <Proj(g), e_i> = <g, e_i> for all e_i in the passed subspace basis.
# Since Proj(g) = sum_j{ a_j e_j } for some a_j's, this gives us a linear system
#   sum_j { <e_j, e_i> a_j } = <g, e_i>
# which we can solve for the a_j.

# Return coefficients for the projection onto the indicated face, relative to the basis elements for the face.
function project_onto_fe_face(g::Function, fe::FENum, face::RelFace, basis::WeakFunsPolyBasis)
  if face == Mesh.interior_face
    const int_mons = WGBasis.interior_mons(basis)
    const ips_bel_vs_bel = WGBasis.ips_interior_mons(fe, basis)
    const ips_g_vs_bels = [Mesh.integral_global_x_face_rel_on_fe_face(g, mon, fe, face, basis.mesh)::R
                           for mon in int_mons]
    ips_bel_vs_bel \ ips_g_vs_bels
  else
    const side_mons = WGBasis.side_mons_for_fe_side(fe, face, basis)
    const ips_bel_vs_bel = WGBasis.ips_fe_side_mons(fe, face, basis)
    const ips_g_vs_bels = [Mesh.integral_global_x_face_rel_on_fe_face(g, mon, fe, face, basis.mesh)::R
                           for mon in side_mons]
    ips_bel_vs_bel \ ips_g_vs_bels
  end
end

function project_onto_fe_face(c::R, fe::FENum, face::RelFace, basis::WeakFunsPolyBasis)
  const mons = face == Mesh.interior_face ? WGBasis.interior_mons(basis) : WGBasis.side_mons_for_fe_side(fe, face, basis)
  const num_mons = length(mons)
  const one_mon = Mesh.one_mon(basis.mesh)
  const proj = zeros(R, num_mons)
  for i=1:num_mons
    if mons[i] == one_mon
      proj[i] = c
      return proj
    end
  end
  error("Could not find one monomial in basis monomials for fe $fe, face $face")
end

project_onto_fe_face_as_poly(g::FunctionOrConst, fe::FENum, face::RelFace, basis::WeakFunsPolyBasis) =
  let proj_coefs = project_onto_fe_face(g, fe, face, basis),
      mons = face == Mesh.interior_face ? WGBasis.interior_mons(basis) : WGBasis.side_mons_for_fe_side(fe, face, basis)
    Polynomial(mons, proj_coefs)
  end


function project_interior_mon_onto_oshape_side(int_mon::Monomial, fe_oshape::OrientedShape, side_face::RelFace, basis::WeakFunsPolyBasis)
  const side_mons = WGBasis.side_mons_for_oshape_side(fe_oshape, side_face, basis)
  const ips_side_bels = WGBasis.ips_oshape_side_mons(fe_oshape, side_face, basis)
  const ips_int_mon_vs_bels = [Mesh.integral_fe_rel_x_side_rel_on_side(int_mon, side_bel_mon, fe_oshape, side_face, basis.mesh)::R
                               for side_bel_mon in side_mons]
  ips_side_bels \ ips_int_mon_vs_bels
end

end # end of module
