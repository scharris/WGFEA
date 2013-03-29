module VBF
export AbstractVariationalBilinearForm,
       is_symmetric,
       interior_bel_vs_interior_bel,
       side_bel_vs_interior_bel,
       interior_bel_vs_side_bel,
       side_bel_vs_side_bel,
       poly_on_fe_face_vs_bel,
       transpose_bel_vs_bel_matrix

import WGBasis.BElNum, WGBasis.WeakFunsPolyBasis
import Poly.Polynomial
import Mesh.FENum, Mesh.FEFace


abstract AbstractVariationalBilinearForm

##########################################################################
# Functions for which methods should be implemented for concrete subtypes.

is_symmetric{BF <: AbstractVariationalBilinearForm}(bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

interior_bel_vs_interior_bel{BF <: AbstractVariationalBilinearForm}(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_bel_vs_interior_bel{BF <: AbstractVariationalBilinearForm}(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

interior_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(bel1::BElNum, bel2::BElNum, basis::WeakFunsPolyBasis, bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")


poly_on_fe_face_vs_bel{BF <: AbstractVariationalBilinearForm}(
  p::Polynomial,
  fe::FENum,
  face::FEFace,
  bel::BElNum,
  basis::WeakFunsPolyBasis,
  bf::BF
) = error("not implemented, bilinear form implementation is incomplete")

#
##########################################################################


function transpose_bel_vs_bel_matrix{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, bf::BF)
  const num_int_bels = basis.num_interior_bels
  const first_side_bel = basis.first_side_bel
  const m = Array(R, basis.total_bels, basis.total_bels)

  if !is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      const ip = interior_bel_vs_interior_bel(bel_num(i), bel_num(j), basis, bf)
      m[j,i] = ip
    end
    for i=first_side_bel:basis.total_bels, j=1:num_interior_bels
      m[j,i] = side_bel_vs_interior_bel(bel_num(i), bel_num(j), basis, bf)
      m[i,j] = side_bel_vs_interior_bel(bel_num(j), bel_num(i), basis, bf)
    end
    for i=first_side_bel:basis.total_bels, j=first_side_bel:basis.total_bels
      m[j,i] = side_bel_vs_side_bel(bel_num(i), bel_num(j), basis, bf)
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_num_i = bel_num(i)
      for j=1:i-1
        const ip = interior_bel_vs_interior_bel(bel_num_i, bel_num(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = interior_bel_vs_interior_bel(bel_num_i, bel_num_i, basis, bf)
    end
    for i=first_side_bel:basis.total_bels, j=1:num_interior_bels
      const ip = side_bel_vs_interior_bel(bel_num(i), bel_num(j), basis, bf)
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
