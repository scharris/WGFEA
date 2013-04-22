module VBF_a_s
export A_s, a_s

using Common
import Poly.Polynomial, Poly.Monomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face, Mesh.OrientedShape
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
                            cache::Dict{Any,Any},
                            bf::A_s)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_interior_mon(monn_1, fe_oshape, basis)
        wgrad_2 = WGBasis.wgrad_interior_mon(monn_2, fe_oshape, basis)
      Mesh.integral_face_rel_on_face(dot(wgrad_1,wgrad_2), fe_oshape, Mesh.interior_face, mesh)
    end

  # stabilization term: The values on the boundary (b_i/j)_b are 0, leaving only
  #   (1/h_T) <Q_b (b_i)_0T, Q_b (b_j)_0T>_bnd(T)
  # = (1/h_T) sum_{s in sides(T)} {<Q_b (b_i)_0T, Q_b (b_j)_0T>_s}
  const stab = begin
    const int_mons = WGBasis.interior_mons(basis)
    const mon_1 = int_mons[monn_1]
    const mon_2 = int_mons[monn_2]
    sum_over_sides = zeroR
    for sf=fe_face(1):Mesh.num_side_faces_for_shape(fe_oshape, mesh)
      const proj_1 = project_interior_mon_onto_oshape_side(mon_1, fe_oshape, sf, basis, cache)
      const proj_2 = project_interior_mon_onto_oshape_side(mon_2, fe_oshape, sf, basis, cache)
      sum_over_sides += Mesh.integral_face_rel_x_face_rel_on_face(proj_1, proj_2, fe_oshape, sf, mesh)
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
                             cache::Dict{Any,Any},
                             bf::A_s)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # weak gradients term
  const ip_wgrads =
    let side_wgrad = WGBasis.wgrad_side_mon(side_monn, fe_oshape, side_face, basis),
        int_wgrad =  WGBasis.wgrad_interior_mon(int_monn, fe_oshape, basis)
      Mesh.integral_face_rel_on_face(dot(side_wgrad, int_wgrad), fe_oshape, Mesh.interior_face, mesh)
    end

  # stabilization term: (1/h_T) <Q_b v_0T - v_b, Q_b w_0T - w_b>_bnd(T)
  # With v being our interior supported basis element and w the side supported element, this becomes
  #    (1/h_T) <Q_b v_0T, -w_b>_bnd(T)
  #  = (1/h_T) <Q_b v_0T, -w_b>_s where s is the support face of w
  const stab = begin
    const side_mon = WGBasis.side_mons_for_fe_side(fe, side_face, basis)[side_monn]
    const int_mon = WGBasis.interior_mons(basis)[int_monn]
    const int_proj = project_interior_mon_onto_oshape_side(int_mon, fe_oshape, side_face, basis, cache)
    const ip = -Mesh.integral_face_rel_x_face_rel_on_face(int_proj, side_mon, fe_oshape, side_face, mesh)
    Mesh.fe_diameter_inv(fe, mesh) * ip
  end

  ip_wgrads + stab
end

import VBF.int_mon_vs_side_mon
function int_mon_vs_side_mon(fe::FENum,
                             int_monn::MonNum,
                             side_monn::MonNum, side_face::FERelFace,
                             basis::WeakFunsPolyBasis,
                             cache::Dict{Any,Any},
                             bf::A_s)
  assert(is_symmetric(bf), "int_mon_vs_side_mon needs independent implementation, bilinear form is not symmetric")
  side_mon_vs_int_mon(fe, side_monn, side_face, int_monn, basis, cache, bf)
end


import VBF.side_mon_vs_side_mon
function side_mon_vs_side_mon(fe::FENum,
                              monn_1::MonNum, side_face_1::FERelFace,
                              monn_2::MonNum, side_face_2::FERelFace,
                              basis::WeakFunsPolyBasis,
                              cache::Dict{Any,Any},
                              bf::A_s)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)

  # weak gradients term
  const ip_wgrads =
    let wgrad_1 = WGBasis.wgrad_side_mon(monn_1, fe_oshape, side_face_1, basis)
        wgrad_2 = WGBasis.wgrad_side_mon(monn_2, fe_oshape, side_face_2, basis)
      Mesh.integral_face_rel_on_face(dot(wgrad_1, wgrad_2), fe_oshape, Mesh.interior_face, mesh)
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
      const side_mons = WGBasis.side_mons_for_fe_side(fe, common_supp_side, basis)
      const mon_1 = side_mons[monn_1]
      const mon_2 = side_mons[monn_2]
      const ip = Mesh.integral_face_rel_x_face_rel_on_face(mon_1, mon_2, fe_oshape, common_supp_side, mesh)
      Mesh.fe_diameter_inv(fe, mesh) * ip
    end

  ip_wgrads + stab
end


function project_interior_mon_onto_oshape_side(mon::Monomial, fe_oshape::OrientedShape, side_face::FERelFace, basis::WeakFunsPolyBasis, cache::Dict{Any,Any})
  const cache_key = (mon, fe_oshape, side_face, basis)
  const cached_proj = get(cache, cache_key, 0)
  if cached_proj != 0
    cached_proj
  else
    const proj = Proj.project_interior_mon_onto_oshape_side(mon, fe_oshape, side_face, basis)
    cache[cache_key] = proj
    proj
  end
end

end # end of module
