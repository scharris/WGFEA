module VBF_a
export A, a

using Common
import Poly.Polynomial, Poly.Monomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face, Mesh.OrientedShape
import Proj
import WGBasis, WGBasis.WeakFunsPolyBasis, WGBasis.BElNum, WGBasis.MonNum
import VBF, VBF.AbstractVariationalBilinearForm

# a(v, w) = sum_T{ (wgrad_T v, wgrad_T w)_T }

type A <: AbstractVariationalBilinearForm
  # projections of interior basis monomials onto fe side faces
  int_mon_side_projs::Array{Array{Array{Polynomial,1},1},1} # by interior monomial number, fe shape, side face

  function A(basis::WeakFunsPolyBasis)
    new(VBF.make_interior_mon_side_projs(basis))
  end
end

a(basis::WeakFunsPolyBasis) = A(basis)


import VBF.is_symmetric
is_symmetric(bf::A) = true


import VBF.int_mon_vs_int_mon
function int_mon_vs_int_mon(fe::FENum,
                            monn_1::MonNum,
                            monn_2::MonNum,
                            basis::WeakFunsPolyBasis,
                            bf::A)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
  const wgrad_1 = WGBasis.wgrad_interior_mon(monn_1, fe_oshape, basis)
  const wgrad_2 = WGBasis.wgrad_interior_mon(monn_2, fe_oshape, basis)
  Mesh.integral_face_rel_on_face(dot(wgrad_1,wgrad_2), fe_oshape, Mesh.interior_face, mesh)
end

import VBF.side_mon_vs_int_mon
function side_mon_vs_int_mon(fe::FENum,
                             side_monn::MonNum, side_face::FERelFace,
                             int_monn::MonNum,
                             basis::WeakFunsPolyBasis,
                             bf::A)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
  const side_wgrad = WGBasis.wgrad_side_mon(side_monn, fe_oshape, side_face, basis)
  const int_wgrad =  WGBasis.wgrad_interior_mon(int_monn, fe_oshape, basis)
  Mesh.integral_face_rel_on_face(dot(side_wgrad, int_wgrad), fe_oshape, Mesh.interior_face, mesh)
end

import VBF.int_mon_vs_side_mon
function int_mon_vs_side_mon(fe::FENum,
                             int_monn::MonNum,
                             side_monn::MonNum, side_face::FERelFace,
                             basis::WeakFunsPolyBasis,
                             bf::A)
  assert(is_symmetric(bf), "int_mon_vs_side_mon needs independent implementation, bilinear form is not symmetric")
  side_mon_vs_int_mon(fe, side_monn, side_face, int_monn, basis, bf)
end


import VBF.side_mon_vs_side_mon
function side_mon_vs_side_mon(fe::FENum,
                              monn_1::MonNum, side_face_1::FERelFace,
                              monn_2::MonNum, side_face_2::FERelFace,
                              basis::WeakFunsPolyBasis,
                              bf::A)
  const mesh = basis.mesh
  const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
  const wgrad_1 = WGBasis.wgrad_side_mon(monn_1, fe_oshape, side_face_1, basis)
  const wgrad_2 = WGBasis.wgrad_side_mon(monn_2, fe_oshape, side_face_2, basis)
  Mesh.integral_face_rel_on_face(dot(wgrad_1, wgrad_2), fe_oshape, Mesh.interior_face, mesh)
end


end # end of module
