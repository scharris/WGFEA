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

using Common
import WGBasis, WGBasis.BElNum, WGBasis.bel_num, WGBasis.WeakFunsPolyBasis, WGBasis.MonNum
import Mesh, Mesh.FENum, Mesh.FERelFace

abstract AbstractVariationalBilinearForm

# We assume the following property for the bilinear forms bf
# represented by subtypes of AbstractVariationalBilinearForm.
# --------------------------------------------------------------------
# (Element Summability)
#   There exist bilinear functions bf_T whose input functions are
#   restricted to domain T, such that
#     bf(v, w) = sum_T { bf_T(v|T, w|T) } for all weak functions v, w
#   where T ranges over all finite elements in the sum.
# --------------------------------------------------------------------
#
# Note that if v or w has no support on a finite elmenent T, then by the
# linearity of bf_T, it follows that bf_T(v|T, w|T) = 0. Thus we can always
# restrict the sum above to just those elements T which include the supports of
# both input functions [Common Support Summability Property]. Further, if v and
# w are weak functions and no finite element includes both supports, then we
# must have bf(v, w) = 0 [Locality Property].
#
# Thus we can always limit our attention to finite elements on which both
# inputs to our bilinear form are non-zero somewhere (not necessarily at any
# common point). And from the summability property, we can compute
# contributions from these finite elements independently and sum them to obtain
# the value for bf. These element-wise functions bf_T are provided by the core
# set of functions [int/side]_mon_vs_[int/side]_mon listed below, which
# together with the is_symmetric function are necessary and sufficient to
# implement variational forms.


##########################################################################
# Functions which should be implemented by specific variational forms.

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
                                                          cache::Dict{Any,Any},
                                                          bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_int_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           int_monn::MonNum,
                                                           basis::WeakFunsPolyBasis,
                                                           cache::Dict{Any,Any},
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

int_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                           int_monn::MonNum,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           basis::WeakFunsPolyBasis,
                                                           cache::Dict{Any,Any},
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe::FENum,
                                                            monn_1::MonNum, side_face_1::FERelFace,
                                                            monn_2::MonNum, side_face_2::FERelFace,
                                                            basis::WeakFunsPolyBasis,
                                                            cache::Dict{Any,Any},
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
# bilinear form bf satisfying the Element Summability requirement.

function int_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, cache::Dict{Any,Any}, bf::BF)
  const fe1 = WGBasis.support_interior_num(bel_1, basis)
  const fe2 = WGBasis.support_interior_num(bel_2, basis)
  if fe1 != fe2
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let monn_1 = WGBasis.interior_mon_num(bel_1, basis)
        monn_2 = WGBasis.interior_mon_num(bel_2, basis)
      int_mon_vs_int_mon(fe1, monn_1, monn_2, basis, cache, bf)
    end
  end
end

function int_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(ibel::BElNum, sbel::BElNum, basis::WeakFunsPolyBasis, cache::Dict{Any,Any}, bf::BF)
  # Determine if the side is one of the faces of the interior's finite element, and which one if so.
  const int_num = WGBasis.support_interior_num(ibel, basis)
  const side_incls = WGBasis.fe_inclusions_of_side_support(sbel, basis)
  const side_face = int_num == side_incls.fe1 ? side_incls.face_in_fe1 :
                    int_num == side_incls.fe2 ? side_incls.face_in_fe2 : Mesh.no_face
  if side_face == Mesh.no_face # no fe includes both supports
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let side_monn = WGBasis.side_mon_num(sbel, basis)
        int_monn = WGBasis.interior_mon_num(ibel, basis)
      int_mon_vs_side_mon(int_num, int_monn, side_monn, side_face, basis, cache, bf)
    end
  end
end


function side_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(sbel::BElNum, ibel::BElNum, basis::WeakFunsPolyBasis, cache::Dict{Any,Any}, bf::BF)
  # Determine if the side is one of the faces of the interior's finite element, and which one if so.
  const side_incls = WGBasis.fe_inclusions_of_side_support(sbel, basis)
  const int_num = WGBasis.support_interior_num(ibel, basis)
  const side_face = int_num == side_incls.fe1 ? side_incls.face_in_fe1 :
                    int_num == side_incls.fe2 ? side_incls.face_in_fe2 : Mesh.no_face
  if side_face == Mesh.no_face
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let side_monn = WGBasis.side_mon_num(sbel, basis)
        int_monn = WGBasis.interior_mon_num(ibel, basis)
      side_mon_vs_int_mon(int_num, side_monn, side_face, int_monn, basis, cache, bf)
    end
  end
end

function side_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, cache::Dict{Any,Any}, bf::BF)
  const incls_1 = WGBasis.fe_inclusions_of_side_support(bel_1, basis)
  const incls_2 = WGBasis.fe_inclusions_of_side_support(bel_2, basis)
  if incls_1.fe1 != incls_2.fe1 &&
     incls_1.fe1 != incls_2.fe2 &&
     incls_1.fe2 != incls_2.fe1 &&
     incls_1.fe2 != incls_2.fe2    # no fe includes both supports...
    zeroR # by locality property
  else
    const monn_1 = WGBasis.side_mon_num(bel_1, basis)
    const monn_2 = WGBasis.side_mon_num(bel_2, basis)

    # By the common support summability property, we can sum the contributions from
    # whichever of the support side including finite elements includes both supports.

    # contribution from incls_1.fe1
    if incls_1.fe1 == incls_2.fe1
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe1, basis, cache, bf)
    elseif incls_1.fe1 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe2, basis, cache, bf)
    else
      zeroR
    end +
    # contribution from incls_1.fe2
    if incls_1.fe2 == incls_2.fe1
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe1, basis, cache, bf)
    elseif incls_1.fe2 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe2, basis, cache, bf)
    else
      zeroR
    end
  end
end



function bel_vs_bel_transpose{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, bf::BF)
  const num_int_bels = basis.num_interior_bels
  const first_nb_side_bel = basis.first_nb_side_bel
  const m = Array(R, basis.total_bels, basis.total_bels)

  const cache = Dict{Any,Any}()

  if !is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      const ip = int_bel_vs_int_bel(bel_num(i), bel_num(j), basis, cache, bf)
      m[j,i] = ip
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const bel_num_i = bel_num(i)
      const bel_num_j = bel_num(j)
      m[j,i] = side_bel_vs_int_bel(bel_num_i, bel_num_j, basis, cache, bf)
      m[i,j] = int_bel_vs_side_bel(bel_num_j, bel_num_i, basis, cache, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=first_nb_side_bel:basis.total_bels
      m[j,i] = side_bel_vs_side_bel(bel_num(i), bel_num(j), basis, cache, bf)
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_num_i = bel_num(i)
      for j=1:i-1
        const ip = int_bel_vs_int_bel(bel_num_i, bel_num(j), basis, cache, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = int_bel_vs_int_bel(bel_num_i, bel_num_i, basis, cache, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const ip = side_bel_vs_int_bel(bel_num(i), bel_num(j), basis, cache, bf)
      m[j,i] = ip
      m[i,j] = ip
    end
    for i=first_nb_side_bel:basis.total_bels
      const bel_num_i = bel_num(i)
      for j=first_nb_side_bel:i-1
        const ip = side_bel_vs_side_bel(bel_num_i, bel_num(j), basis, cache, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = side_bel_vs_side_bel(bel_num_i, bel_num_i, basis, cache, bf)
    end
  end
  m
end

end # end of module
