module VBF_a_s
export VBF

using Common
import Mesh, Mesh.FENum, Mesh.FEFace, Mesh.fe_face
import Poly.Polynomial
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.mon_num, WGBasis.bel_num
import WGrad
import VBF.AbstractVariationalBilinearForm


type VBF <: AbstractVariationalBilinearForm
end


# Change to false if a post-multiplication matrix ("a") is introduced for the left wgrad in the inner product
is_symmetric{BF <: AbstractVariationalBilinearForm}(bf::BF) = true

function interior_bel_vs_interior_bel(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::VBF)
  const fe1 = WGBasis.support_interior_num(bel1, basis)
  const fe2 = WGBasis.support_interior_num(bel2, basis)
  if fe1 != fe2
    zeroR
  else
    const mesh = basis.mesh
    const mon_num_1 = WGBasis.interior_monomial_num(bel1, basis)
    const mon_num_2 = WGBasis.interior_monomial_num(bel2, basis)

    # weak gradients term
    const ip_wgrads =
      let wgrad_1 = WGBasis.wgrad_for_interior_mon_num(mon_num_1, basis) # TODO: Would post-multiply by matrix "a" here.
          wgrad_2 = WGBasis.wgrad_for_interior_mon_num(mon_num_2, basis)
        Mesh.integral_on_ref_fe_face(dot(wgrad_1, wgrad_2), Mesh.interior_face, mesh)
      end

    # stabilization term: The values on the actual boundary (b_i/j)_b are 0, leaving: <b_i, b_j>_boundary(T)
    const stab = begin
      const mon_1 = WGBasis.interior_monomial_by_num(mon_num_1, basis)
      const mon_2 = WGBasis.interior_monomial_by_num(mon_num_2, basis)
      sum_over_sides = zeroR
      for s=1:Mesh.num_side_faces_per_fe(mesh)
        const side = fe_face(s)
        sum_over_sides += Mesh.integral_prod_on_ref_fe_face(mon_1, mon_2, side, mesh)
      end
      sum_over_sides
    end

    ip_wgrads + stab
  end
end

function interior_bel_vs_side_bel(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::VBF)
# TODO
end

function side_bel_vs_interior_bel(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::VBF)
# TODO
end

function side_bel_vs_side_bel(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::VBF)
# TODO
end

poly_on_fe_face_vs_bel(
  p::Polynomial,
  fe::FENum,
  face::FEFace,
  bel::BElNum,
  basis::WeakFunsPolyBasis,
  bf::VBF
) = "TODO"



end # end of module
