module VBF_a_s
export A_s, a_s

using Common
import Poly.Polynomial
import Mesh, Mesh.OShapeNum, Mesh.FEFaceNum, Mesh.feface_one
import Proj
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.MonNum
import VBF, VBF.AbstractVariationalBilinearForm

# a_s(v, w) = sum_T{ (a wgrad_T v, wgrad_T w)_T } + sum_T{ (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T) }
# For now we will let the post-multiplying matrix a in the first term be the identity matrix.

immutable A_s <: AbstractVariationalBilinearForm
  # projections of interior basis monomials onto fe side faces
  int_mon_side_projs::Array{Array{Array{Polynomial,1},1},1} # by interior monomial number, fe shape, side face

  function A_s(basis::WeakFunsPolyBasis)
    new(Proj.make_interior_mon_side_projs(basis))
  end
end

a_s(basis::WeakFunsPolyBasis) = A_s(basis)


# Change to false if a non-identity post-multiplying matrix "a" is introduced for the left wgrad in the inner product
import VBF.is_symmetric
is_symmetric(bf::A_s) = true


import VBF.int_mon_vs_int_mon
function int_mon_vs_int_mon(fe_oshape::OShapeNum,
                            monn_1::MonNum,
                            monn_2::MonNum,
                            basis::WeakFunsPolyBasis,
                            bf::A_s)
  const mesh = basis.mesh

  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_interior_mon(monn_1, fe_oshape, basis)
        wgrad_2 = WGBasis.wgrad_interior_mon(monn_2, fe_oshape, basis)
      Mesh.integral_face_rel_on_oshape_face(dot(wgrad_1,wgrad_2), fe_oshape, Mesh.interior_face, mesh)
    end

  # stabilization term: The values on the boundary (b_i/j)_b are 0, leaving only
  #   (1/h_T) <Q_b (b_i)_0T, Q_b (b_j)_0T>_bnd(T)
  # = (1/h_T) sum_{s in sides(T)} {<Q_b (b_i)_0T, Q_b (b_j)_0T>_s}
  const stab = begin
    sum_over_sides = zeroR
    for sf=feface_one:Mesh.num_side_faces_for_shape(fe_oshape, mesh)
      const proj_1 = bf.int_mon_side_projs[monn_1][fe_oshape][sf]
      const proj_2 = bf.int_mon_side_projs[monn_2][fe_oshape][sf]
      sum_over_sides += Mesh.integral_face_rel_x_face_rel_on_oshape_face(proj_1, proj_2, fe_oshape, sf, mesh)
    end
    Mesh.shape_diameter_inv(fe_oshape, mesh) * sum_over_sides
  end

  ip_wgrads + stab
end

import VBF.side_mon_vs_int_mon
function side_mon_vs_int_mon(fe_oshape::OShapeNum,
                             side_monn::MonNum, side_face::FEFaceNum,
                             int_monn::MonNum,
                             basis::WeakFunsPolyBasis,
                             bf::A_s)
  const mesh = basis.mesh

  # weak gradients term
  const ip_wgrads =
    let side_wgrad = WGBasis.wgrad_side_mon(side_monn, fe_oshape, side_face, basis),
        int_wgrad =  WGBasis.wgrad_interior_mon(int_monn, fe_oshape, basis)
      Mesh.integral_face_rel_on_oshape_face(dot(side_wgrad, int_wgrad), fe_oshape, Mesh.interior_face, mesh)
    end

  # stabilization term: (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T)
  # With v being our interior supported basis element and w the side supported element, this becomes
  #    (1/h_T) <Q_b v_0T, -w_b>_bnd(T)
  #  = (1/h_T) <Q_b v_0T, -w_b>_s where s is the support face of w
  const stab = begin
    const side_mon = WGBasis.side_mons_for_oshape_side(fe_oshape, side_face, basis)[side_monn]
    const int_proj = bf.int_mon_side_projs[int_monn][fe_oshape][side_face]
    const ip = -Mesh.integral_face_rel_x_face_rel_on_oshape_face(int_proj, side_mon, fe_oshape, side_face, mesh)
    Mesh.shape_diameter_inv(fe_oshape, mesh) * ip
  end

  ip_wgrads + stab
end

import VBF.int_mon_vs_side_mon
function int_mon_vs_side_mon(fe_oshape::OShapeNum,
                             int_monn::MonNum,
                             side_monn::MonNum, side_face::FEFaceNum,
                             basis::WeakFunsPolyBasis,
                             bf::A_s)
  assert(is_symmetric(bf), "int_mon_vs_side_mon needs independent implementation, bilinear form is not symmetric")
  side_mon_vs_int_mon(fe_oshape, side_monn, side_face, int_monn, basis, bf)
end


import VBF.side_mon_vs_side_mon
function side_mon_vs_side_mon(fe_oshape::OShapeNum,
                              monn_1::MonNum, side_face_1::FEFaceNum,
                              monn_2::MonNum, side_face_2::FEFaceNum,
                              basis::WeakFunsPolyBasis,
                              bf::A_s)
  const mesh = basis.mesh

  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_side_mon(monn_1, fe_oshape, side_face_1, basis)
        wgrad_2 = WGBasis.wgrad_side_mon(monn_2, fe_oshape, side_face_2, basis)
      Mesh.integral_face_rel_on_oshape_face(dot(wgrad_1, wgrad_2), fe_oshape, Mesh.interior_face, mesh)
    end

  # stabilization term: (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T)
  # With v and w being our side supported elements, this becomes
  #    (1/h_T) <-v_b, -w_b>_bnd(T)
  #  = (1/h_T) sum_{s in sides(T)} <v_b, w_b>_s
  #  = | (1/h_T) <v_b, w_b>_s,  if v and w have a common support side face s
  #    | 0, otherwise
  const stab =
    if side_face_1 != side_face_2
      zeroR
    else
      const common_supp_side = side_face_1
      const side_mons = WGBasis.side_mons_for_oshape_side(fe_oshape, common_supp_side, basis)
      const mon_1 = side_mons[monn_1]
      const mon_2 = side_mons[monn_2]
      const ip = Mesh.integral_face_rel_x_face_rel_on_oshape_face(mon_1, mon_2, fe_oshape, common_supp_side, mesh)
      Mesh.shape_diameter_inv(fe_oshape, mesh) * ip
    end

  ip_wgrads + stab
end


end # end of module
