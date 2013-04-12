module VBF_a_s
export A_s, a_s

using Common
import Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.MonNum
import Proj
import VBF.AbstractVariationalBilinearForm

# a_s(v, w) = sum_T{ (a wgrad_T v, wgrad_T w)_T } + sum_T{ (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T) }
# For now we will let the post-multiplying matrix a in the first term be the identity matrix.

type A_s <: AbstractVariationalBilinearForm
end

const a_s = A_s()

# Change to false if a non-identity post-multiplying matrix "a" is introduced for the left wgrad in the inner product
import VBF.is_symmetric
is_symmetric(bf::A_s) = true


import VBF.int_mon_vs_int_mon
function int_mon_vs_int_mon(fe::FENum,
                            monn_1::MonNum,
                            monn_2::MonNum,
                            basis::WeakFunsPolyBasis,
                            bf::A_s)
  const mesh = basis.mesh

  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_for_mon_on_interior(monn_1, basis)
        wgrad_2 = WGBasis.wgrad_for_mon_on_interior(monn_2, basis)
      Mesh.integral_face_rel_on_face(dot(wgrad_1,wgrad_2), Mesh.interior_face, mesh)
    end

  # stabilization term: The values on the boundary (b_i/j)_b are 0, leaving only
  #   (1/h_T) <Q_b (b_i)_0T, Q_b (b_j)_0T>_bnd(T)
  # = (1/h_T) sum_{s in sides(T)} {<Q_b (b_i)_0T, Q_b (b_j)_0T>_s}
  const stab = begin
    const mon_1 = WGBasis.int_monomial_by_num(monn_1, basis)
    const mon_2 = WGBasis.int_monomial_by_num(monn_2, basis)
    sum_over_sides = zeroR
    for s=fe_face(1):fe_face(Mesh.num_side_faces_per_fe(mesh))
      const proj_1 = Proj.project_interior_monomial_onto_side_face(mon_1, s, basis)
      const proj_2 = Proj.project_interior_monomial_onto_side_face(mon_2, s, basis)
      sum_over_sides += Mesh.integral_face_rel_x_face_rel_on_face(proj_1, proj_2, s, mesh)
    end
    Mesh.fe_diameter_inv(fe, mesh) * sum_over_sides
  end

  ip_wgrads + stab
end

import VBF.side_mon_vs_int_mon
function side_mon_vs_int_mon(fe::FENum,
                             side_monn::MonNum, side_face::FERelFace,
                             int_monn::MonNum,
                             basis::WeakFunsPolyBasis,
                             bf::A_s)
  # weak gradients term
  const ip_wgrads =
    let side_wgrad = WGBasis.wgrad_for_mon_on_side_face(side_monn, side_face, basis)
        int_wgrad =  WGBasis.wgrad_for_mon_on_interior(int_monn, basis),
      Mesh.integral_face_rel_on_face(dot(side_wgrad, int_wgrad), Mesh.interior_face, mesh)
    end

  # stabilization term: (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T)
  # With v being our interior supported basis element and w the side supported element, this becomes
  #    (1/h_T) <Q_b v_0T, -w_b>_bnd(T)
  #  = (1/h_T) (-<Q_b v_0T, w_b>_s) where s is the support face of w
  const stab = begin
    const side_mon = WGBasis.side_monomial_by_face_and_num(side_face, side_monn, basis)
    const int_mon = WGBasis.int_monomial_by_num(int_monn, basis)
    const int_proj = Proj.project_interior_monomial_onto_side_face(int_mon, side_face, basis)
    const ip = -Mesh.integral_face_rel_x_face_rel_on_face(int_proj, side_mon, side_face, mesh)
    Mesh.fe_diameter_inv(fe, mesh) * ip
  end

  ip_wgrads + stab
end

import VBF.int_mon_vs_side_mon
function int_mon_vs_side_mon(fe::FENum,
                             int_monn::MonNum,
                             side_monn::MonNum, side_face::FERelFace,
                             basis::WeakFunsPolyBasis,
                             bf::A_s)
  assert(is_symmetric(bf), "int_mon_vs_side_mon needs independent implementation, bilinear form is not symmetric")
  side_mon_vs_int_mon(fe, side_monn, side_face, int_monn, basis, bf)
end


import VBF.side_mon_vs_side_mon
function side_mon_vs_side_mon(fe::FENum,
                              monn_1::MonNum, side_face_1::FERelFace,
                              monn_2::MonNum, side_face_2::FERelFace,
                              basis::WeakFunsPolyBasis,
                              bf::A_s)
  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_for_mon_on_side_face(monn_1, side_face_1, basis)
        wgrad_2 = WGBasis.wgrad_for_mon_on_side_face(monn_2, side_face_2, basis)
      Mesh.integral_face_rel_on_face(dot(wgrad_1, wgrad_2), Mesh.interior_face, basis.mesh)
    end

  # stabilization term: (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T)
  # With v and w being our side supported elements, this becomes
  #    (1/h_T) <-v_b, -w_b>_bnd(T)
  #  = (1/h_T) sum_{s in sides(T)} <v_b, w_b>_s
  #  = | (1/h_T) <v_b, w_b>_s',  if v and w have a common support side face s'
  #    | 0, otherwise
  const stab =
    if side_face_1 != side_face_2
      zeroR
    else
      const common_supp_side = side_face_1
      const mon_1 = WGBasis.side_monomial_by_face_and_num(side_face_1, monn_1, basis)
      const mon_2 = WGBasis.side_monomial_by_face_and_num(side_face_2, monn_2, basis)
      const ip = Mesh.integral_face_rel_x_face_rel_on_face(mon_1, mon_2, common_supp_side, mesh)
      Mesh.fe_diameter_inv(fe, mesh) * ip
    end

  ip_wgrads + stab
end

end # end of module
