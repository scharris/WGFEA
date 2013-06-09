module Proj
export project_onto_fe_face_supported_approx_subspace,
       project_onto_fe_face_supported_approx_subspace_as_poly,
       project_interior_mon_onto_oshape_side_supported_approx_subspace

using Common
import Poly.Monomial, Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.feface_one, Mesh.oshape_one
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.monnum


# Projection onto a Subspace
# --------------------------
# <Proj(g), e_i> = <g, e_i> for all e_i in the passed subspace basis.
# Since Proj(g) = sum_j{ a_j e_j } for some a_j's, this gives us a linear system
#   sum_j { <e_j, e_i> a_j } = <g, e_i>
# which we can solve for the a_j.


# Return coefficients for the projection onto the indicated face, relative to the basis elements for the face.
function project_onto_fe_face_supported_approx_subspace(g::Function,
                                                        fe::FENum, face::FEFaceNum,
                                                        basis::WeakFunsPolyBasis)
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


function project_onto_fe_face_supported_approx_subspace(c::R,
                                                        fe::FENum, face::FEFaceNum,
                                                        basis::WeakFunsPolyBasis)
  const mons = face == Mesh.interior_face ? WGBasis.interior_mons(basis) :
                                            WGBasis.side_mons_for_fe_side(fe, face, basis)
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



function project_onto_fe_face_supported_approx_subspace_as_poly(g::FunctionOrConst,
                                                                fe::FENum, face::FEFaceNum,
                                                                basis::WeakFunsPolyBasis)
  const proj_coefs = project_onto_fe_face_supported_approx_subspace(g, fe, face, basis)
  const mons = face == Mesh.interior_face ? WGBasis.interior_mons(basis) :
                                            WGBasis.side_mons_for_fe_side(fe, face, basis)
  Polynomial(mons, proj_coefs)
end

# TODO: unit tests
function make_interior_mon_side_projs(basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const int_mons = WGBasis.interior_mons(basis)
  const num_int_mons = monnum(length(int_mons))
  const projs_by_int_monn = Array(Array{Array{Polynomial,1},1}, num_int_mons)
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  for int_monn=monnum(1):num_int_mons
    projs_by_oshape = Array(Array{Polynomial,1}, num_oshapes)
    for os=oshape_one:num_oshapes
      const sides_per_fe = Mesh.num_side_faces_for_shape(os, mesh)
      const projs_by_side = Array(Polynomial, sides_per_fe)
      for sf=feface_one:sides_per_fe
        const side_mons = WGBasis.side_mons_for_oshape_side(os, sf, basis)
        const proj_coefs = project_interior_mon_onto_oshape_side_supported_approx_subspace(int_mons[int_monn], os, sf, basis)
        projs_by_side[sf] = Polynomial(side_mons, proj_coefs)
      end
      projs_by_oshape[os] = projs_by_side
    end
    projs_by_int_monn[int_monn] = projs_by_oshape
  end
  projs_by_int_monn
end

function project_interior_mon_onto_oshape_side_supported_approx_subspace(int_mon::Monomial,
                                                                         fe_oshape::OShapeNum, side_face::FEFaceNum,
                                                                         basis::WeakFunsPolyBasis)
  const side_mons = WGBasis.side_mons_for_oshape_side(fe_oshape, side_face, basis)
  const ips_side_bels = WGBasis.ips_oshape_side_mons(fe_oshape, side_face, basis)
  const ips_int_mon_vs_bels =
    [Mesh.integral_fe_rel_x_side_rel_on_oshape_side(int_mon, side_bel_mon, fe_oshape, side_face, basis.mesh)::R
     for side_bel_mon in side_mons]

  ips_side_bels \ ips_int_mon_vs_bels
end

end # end of module
