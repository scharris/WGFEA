module VBF
export AbstractVariationalBilinearForm,
       is_symmetric,
       int_mon_vs_int_mon,
       side_mon_vs_int_mon,
       int_mon_vs_side_mon,
       side_mon_vs_side_mon,
       int_bel_vs_int_bel,
       side_bel_vs_int_bel,
       int_bel_vs_side_bel,
       side_bel_vs_side_bel,
       bel_vs_bel_transpose

import WGBasis.BElNum, WGBasis.WeakFunsPolyBasis, WGBasis.MonNum
import Mesh.FENum, Mesh.FERelFace


abstract AbstractVariationalBilinearForm

##########################################################################
# Functions which should be implemented as methods for concrete subtypes.

is_symmetric{BF <: AbstractVariationalBilinearForm}(bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

# bilinear form evaluation functions

# In all of the following evaluation functions, monomials are interpreted as
# having face-relative origins as determined by the mesh embedded in the basis.
# That is, they are implictly precomposed with a translation x -> x - x_0(F)
# where x_0(F) depends on the interior or side in a way determined by the mesh.
# Implementations which use mesh integration functions should make sure that
# the functions called interpret monomials in the expected way.

int_mon_vs_int_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                          monn_1::MonNum,
                                                          monn_2::MonNum,
                                                          basis::WeakFunsPolyBasis,
                                                          bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_int_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           int_monn::MonNum,
                                                           basis::WeakFunsPolyBasis,
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

int_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                           int_monn::MonNum,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           basis::WeakFunsPolyBasis,
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                            monn_1::MonNum, side_face_1::FERelFace,
                                                            monn_2::MonNum, side_face_2::FERelFace,
                                                            basis::WeakFunsPolyBasis,
                                                            bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

#
##########################################################################


# Functions with Default Implementations
#
# The functions below should not need to be implemented by subtypes, unless
# to take advantage of special shortcut evaluations for efficiency that
# may be possible for the subtype. These default methods rely on the required
# functions above. The implementations should be valid for any variational
# bilinear form bf satisfying the following locality assumption:
# (Locality Assumption)
#   For every pair of weak functions v and w, if there is no finite element
#   which includes the supports of both v and w, then bf(v, w) = 0.


function int_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  const fe1 = WGBasis.support_interior_num(bel_1, basis)
  const fe2 = WGBasis.support_interior_num(bel_2, basis)
  if fe1 != fe2
    zeroR
  else
    let monn_1 = WGBasis.int_monomial_num(bel_1, basis)
        monn_2 = WGBasis.int_monomial_num(bel_2, basis)
      int_mon_vs_int_mon(fe1, monn_1, monn_2, basis, bf)
    end
  end
end

function int_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(ibel::BElNum, sbel::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  # Determine if the side is one of the faces of the interior's finite element, and which one if so.
  const int_num = WGBasis.support_interior_num(ibel, basis)
  const side_incls = WGBasis.fe_inclusions_of_side_support(sbel, basis)
  const side_face = int_num == side_incls.fe1 ? side_incls.face_in_fe1 :
                    int_num == side_incls.fe2 ? side_incls.face_in_fe2 : Mesh.no_face
  if side_face == Mesh.no_face
    zeroR
  else
    let side_mon_num = WGBasis.side_monomial_num(sbel, basis)
        int_monn = WGBasis.int_monomial_num(ibel, basis)
      int_mon_vs_side_mon(int_num, int_monn, side_mon_num, side_face, basis, bf)
    end
  end
end


function side_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(sbel::BElNum, ibel::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  # Determine if the side is one of the faces of the interior's finite element, and which one if so.
  const side_incls = WGBasis.fe_inclusions_of_side_support(sbel, basis)
  const int_num = WGBasis.support_interior_num(ibel, basis)
  const side_face = int_num == side_incls.fe1 ? side_incls.face_in_fe1 :
                    int_num == side_incls.fe2 ? side_incls.face_in_fe2 : Mesh.no_face
  if side_face == Mesh.no_face
    zeroR
  else
    let side_mon_num = WGBasis.side_monomial_num(sbel, basis)
        int_monn = WGBasis.int_monomial_num(ibel, basis)
      side_mon_vs_int_mon(int_num, side_mon_num, side_face, int_monn, basis, bf)
    end
  end
end

function side_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  const incls_1 = WGBasis.fe_inclusions_of_side_support(bel_1, basis)
  const incls_2 = WGBasis.fe_inclusions_of_side_support(bel_2, basis)
  if incls_1.fe1 != incls_2.fe1 &&
     incls_1.fe1 != incls_2.fe2 &&
     incls_1.fe2 != incls_2.fe1 &&
     incls_1.fe2 != incls_2.fe2    # no fe including both supports
    zeroR
  else
    const monn_1 = WGBasis.side_monomial_num(bel_1, basis)
    const monn_2 = WGBasis.side_monomial_num(bel_2, basis)

    # contribution from incls_1.fe1
    if incls_1.fe1 == incls_2.fe1
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe1, bf)
    elseif incls_1.fe1 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe2, bf)
    else
      zeroR
    end +
    # contribution from incls_1.fe2
    if incls_1.fe2 == incls_2.fe1
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe1, bf)
    elseif incls_1.fe2 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe2, bf)
    else
      zeroR
    end
  end
end



function bel_vs_bel_transpose{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, bf::BF)
  const num_int_bels = basis.num_int_bels
  const first_side_bel = basis.first_side_bel
  const m = Array(R, basis.total_bels, basis.total_bels)

  if !is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      const ip = int_bel_vs_int_bel(bel_num(i), bel_num(j), basis, bf)
      m[j,i] = ip
    end
    for i=first_side_bel:basis.total_bels, j=1:num_int_bels
      const bel_num_i = bel_num(i)
      const bel_num_j = bel_num(j)
      m[j,i] = side_bel_vs_int_bel(bel_num_i, bel_num_j, basis, bf)
      m[i,j] = int_bel_vs_side_bel(bel_num_j, bel_num_i, basis, bf)
    end
    for i=first_side_bel:basis.total_bels, j=first_side_bel:basis.total_bels
      m[j,i] = side_bel_vs_side_bel(bel_num(i), bel_num(j), basis, bf)
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_num_i = bel_num(i)
      for j=1:i-1
        const ip = int_bel_vs_int_bel(bel_num_i, bel_num(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = int_bel_vs_int_bel(bel_num_i, bel_num_i, basis, bf)
    end
    for i=first_side_bel:basis.total_bels, j=1:num_int_bels
      const ip = side_bel_vs_int_bel(bel_num(i), bel_num(j), basis, bf)
      m[j,i] = ip
      m[i,j] = ip
    end
    for i=first_side_bel:basis.total_bels
      const bel_num_i = bel_num(i)
      for j=first_side_bel:i-1
        const ip = side_bel_vs_side_bel(bel_num_i, bel_num(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = side_bel_vs_side_bel(bel_num_i, bel_num_i, basis, bf)
    end
  end
  m
end

end # end of module
