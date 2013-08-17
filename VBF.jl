module VBF
export AbstractVariationalBilinearForm,
       is_symmetric,
       int_mon_vs_int_mon,
       side_mon_vs_int_mon,
       int_mon_vs_side_mon,
       side_mon_vs_side_mon,
       poly_on_face_vs_poly_on_face,
       bel_vs_bel_transpose

require("Common.jl")
require("Poly.jl")
require("Mesh.jl")
require("RMesh.jl")
require("WGrad.jl")
require("WGBasis.jl")
require("Proj.jl")
require("ParCtrl.jl")

using Common
import Poly.Polynomial, Poly.Monomial
import Mesh, Mesh.AbstractMesh, Mesh.OShapeNum, Mesh.FEFaceNum, Mesh.feface_one, Mesh.fefacenum, Mesh.oshape_one, Mesh.fenum
import Proj
import WGBasis, WGBasis.BElNum, WGBasis.WeakFunsPolyBasis, WGBasis.MonNum, WGBasis.monnum
import ParCtrl

abstract AbstractVariationalBilinearForm

# TODO: Rewrite this to document the even stronger assumption being used
#       in the code, which is that all the bf_T's for T's of the same oriented
#       shape are the *same* function after translation of the input functions,
#       ie. they are all implementable as a single bf per reference element.
#
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

is_symmetric(bf::AbstractVariationalBilinearForm) =
  error("not implemented, bilinear form implementation is incomplete")

# bilinear form evaluation functions

# In all of the following evaluation functions, monomials are interpreted as
# having face-relative origins as determined by the mesh embedded in the basis.
# That is, they are implictly precomposed with a translation x -> x - x_0(F)
# where x_0(F) depends on the interior or side in a way determined by the mesh.
# Implementations which use mesh integration functions should make sure that
# the functions called interpret monomials in the expected way.

int_mon_vs_int_mon(fe_oshape::OShapeNum,
                   monn_1::MonNum,
                   monn_2::MonNum,
                   basis::WeakFunsPolyBasis,
                   bf::AbstractVariationalBilinearForm) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_int_mon(fe_oshape::OShapeNum,
                    side_monn::MonNum, side_face::FEFaceNum,
                    int_monn::MonNum,
                    basis::WeakFunsPolyBasis,
                    bf::AbstractVariationalBilinearForm) =
  error("not implemented, bilinear form implementation is incomplete")

int_mon_vs_side_mon(fe_oshape::OShapeNum,
                    int_monn::MonNum,
                    side_monn::MonNum, side_face::FEFaceNum,
                    basis::WeakFunsPolyBasis,
                    bf::AbstractVariationalBilinearForm) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_side_mon(fe_oshape::OShapeNum,
                     monn_1::MonNum, side_face_1::FEFaceNum,
                     monn_2::MonNum, side_face_2::FEFaceNum,
                     basis::WeakFunsPolyBasis,
                     bf::AbstractVariationalBilinearForm) =
  error("not implemented, bilinear form implementation is incomplete")

#
##########################################################################


function poly_on_face_vs_poly_on_face(fe_oshape::OShapeNum,
                                      p1_coefs::Array{R,1}, face_1::FEFaceNum,
                                      p2_coefs::Array{R,1}, face_2::FEFaceNum,
                                      basis::WeakFunsPolyBasis,
                                      bf::AbstractVariationalBilinearForm)
  term_pairs_sum = zeroR
  if face_1 == Mesh.interior_face # interior vs interior or side
    if face_2 == Mesh.interior_face # interior vs interior
      const num_int_mons = WGBasis.mons_per_fe_interior(basis)
      assert(length(p1_coefs) == length(p2_coefs) == num_int_mons)
      for monn_1=monnum(1):num_int_mons, monn_2=monnum(1):num_int_mons
        term_pairs_sum += p1_coefs[monn_1]*p2_coefs[monn_2] *
                          int_mon_vs_int_mon(fe_oshape, monn_1, monn_2, basis, bf)
      end
    else # interior vs side
      const num_int_mons, num_side_mons = WGBasis.mons_per_fe_interior(basis), WGBasis.mons_per_fe_side(basis)
      assert(length(p1_coefs) == num_int_mons && length(p2_coefs) == num_side_mons)
      for monn_1=monnum(1):num_int_mons, monn_2=monnum(1):num_side_mons
        term_pairs_sum += p1_coefs[monn_1]*p2_coefs[monn_2] *
                          int_mon_vs_side_mon(fe_oshape, monn_1, monn_2, face_2, basis, bf)
      end
    end
  else # side vs interior or side
    if face_2 == Mesh.interior_face # side vs interior
      const num_int_mons, num_side_mons = WGBasis.mons_per_fe_interior(basis), WGBasis.mons_per_fe_side(basis)
      assert(length(p1_coefs) == num_side_mons && length(p2_coefs) == num_int_mons)
      for monn_1=monnum(1):num_side_mons, monn_2=monnum(1):num_int_mons
        term_pairs_sum += p1_coefs[monn_1]*p2_coefs[monn_2] *
                          side_mon_vs_int_mon(fe_oshape, monn_1, face_1, monn_2, basis, bf)
      end
    else # side vs side
      const num_side_mons = WGBasis.mons_per_fe_side(basis)
      assert(length(p1_coefs) == length(p2_coefs) == num_side_mons)
      for monn_1=monnum(1):num_side_mons, monn_2=monnum(1):num_side_mons
        term_pairs_sum += p1_coefs[monn_1]*p2_coefs[monn_2] *
                          side_mon_vs_side_mon(fe_oshape, monn_1, face_1, monn_2, face_2, basis, bf)
      end
    end
  end
  term_pairs_sum
end


function bel_vs_bel_transpose(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
  # data arrays for construction of the sparse matrix
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  const row_nums = Array(Int, interacting_bel_pairs_ub)
  const col_nums = Array(Int, interacting_bel_pairs_ub)
  const nonzeros = Array(R, interacting_bel_pairs_ub)

  const mesh = basis.mesh
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)

  # Precompute vbf values for each oriented shape and pairing of monomial/support face with monomial/support face.
  const int_int_vbf_vals, side_int_vbf_vals, side_side_vbf_vals, int_side_vbf_vals =
    if ParCtrl.parallel_basis_vbf_vs_vbf_vals
        const vals = pmap(f -> f(basis, vbf),
                          (ref_int_vs_int_vbf_values, ref_side_vs_int_vbf_vals, ref_side_vs_side_vbf_vals))
        const int_vs_side = is_symmetric(vbf) ? vals[2] : ref_int_vs_side_vbf_vals(basis, vbf)
        vals[1], vals[2], vals[3], int_vs_side
    else
      const side_vs_int = ref_side_vs_int_vbf_vals(basis, vbf)
      (ref_int_vs_int_vbf_values(basis, vbf), 
       side_vs_int,
       ref_side_vs_side_vbf_vals(basis, vbf),
       is_symmetric(vbf) ? side_vs_int : ref_int_vs_side_vbf_vals(basis, vbf))
    end

  # work array for remembering which sides are nb sides within a finite element
  const is_nb_side = Array(Bool, Mesh.max_num_shape_sides(mesh))

  nnz = 0 # current number of non-zero values stored, the last stored position in the data arrays

  for fe=fenum(1):Mesh.num_fes(mesh)
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
    const fe_num_sides = Mesh.num_side_faces_for_shape(fe_oshape, mesh)

    # fill is_nb_side work array (up to this shape's # sides)
    for sf=feface_one:fe_num_sides  is_nb_side[sf] = !Mesh.is_boundary_side(fe, sf, mesh) end

    # fill interior vs interior matrix values
    for monn_1=monnum(1):num_int_mons, monn_2=monnum(1):num_int_mons
      const vbf_val = int_int_vbf_vals[fe_oshape][monn_1, monn_2]
      if vbf_val != zeroR
        const beln_1 = WGBasis.interior_mon_bel_num(fe, monn_1, basis)
        const beln_2 = WGBasis.interior_mon_bel_num(fe, monn_2, basis)
        nnz += 1
        nonzeros[nnz] = vbf_val
        row_nums[nnz] = beln_2 # 2nd first because result matrix is transpose of bel vs bel matrix
        col_nums[nnz] = beln_1
      end
    end

    # fill remaining matrix values, which all involve at least one side basis element
    for sf=feface_one:fe_num_sides if is_nb_side[sf]
      for sf_monn=monnum(1):num_side_mons
        const sf_beln = WGBasis.side_mon_bel_num(fe, sf, sf_monn, basis)

        # fill side vs interior and interior vs side matrix values
        for int_monn=monnum(1):num_int_mons
          int_beln = nothing::Union(BElNum, Nothing)
          # side vs interior
          const side_int_vbf_val = side_int_vbf_vals[fe_oshape][sf, sf_monn, int_monn]
          if side_int_vbf_val != zeroR
            int_beln = WGBasis.interior_mon_bel_num(fe, int_monn, basis)
            nnz += 1
            nonzeros[nnz] = side_int_vbf_val
            row_nums[nnz] = int_beln # 2nd first because result matrix is transpose of bel vs bel matrix
            col_nums[nnz] = sf_beln
          end
          # interior vs side
          const int_side_vbf_val = int_side_vbf_vals[fe_oshape][sf, sf_monn, int_monn]
          if int_side_vbf_val != zeroR
            if int_beln == nothing
              int_beln = WGBasis.interior_mon_bel_num(fe, int_monn, basis)
            end
            nnz += 1
            nonzeros[nnz] = int_side_vbf_val
            row_nums[nnz] = sf_beln # 2nd first because result matrix is transpose of bel vs bel matrix
            col_nums[nnz] = int_beln
          end
        end # interiors (with sf)

        # fill side vs side matrix elements
        for sf2=feface_one:fe_num_sides if is_nb_side[sf2]
          for sf2_monn=monnum(1):num_side_mons
            const vbf_val = side_side_vbf_vals[fe_oshape][sf, sf_monn, sf2, sf2_monn]
            if vbf_val != zeroR
              const sf2_beln = WGBasis.side_mon_bel_num(fe, sf2, sf2_monn, basis)
              nnz += 1
              nonzeros[nnz] = vbf_val
              row_nums[nnz] = sf2_beln # 2nd first because result matrix is transpose of bel vs bel matrix
              col_nums[nnz] = sf_beln
            end
          end # sf2 mons
        end end # second nb side faces, sf2

      end # first side mons
    end end # first nb side faces
  end # fes

  sparse(row_nums[1:nnz],
         col_nums[1:nnz],
         nonzeros[1:nnz],
         basis.total_bels,
         basis.total_bels)
end

function ref_int_vs_int_vbf_values(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
  const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
  const int_int_vbf_vals = Array(Array{R,2}, num_oshapes) # vbf values indexed by fe oshape, (monn_1, monn_2)
  for os=oshape_one:num_oshapes
    int_int_vbf_vals[os] = int_vs_int_vbf_vals_for_oshape(os, basis, vbf)
  end
  int_int_vbf_vals
end

function ref_side_vs_int_vbf_vals(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
  const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
  const side_int_vbf_vals  = Array(Array{R,3}, num_oshapes) # by fe oshape, (side_face, side_monn, int_monn)
  for os=oshape_one:num_oshapes
    side_int_vbf_vals[os] = side_vs_int_vbf_vals_for_oshape(os, basis, vbf)
  end
  side_int_vbf_vals
end

function ref_int_vs_side_vbf_vals(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
  const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
  const int_side_vbf_vals  = Array(Array{R,3}, num_oshapes) # by fe oshape, (side_face, side_monn, int_monn) # (ordered for compat with side-vs-int array)
  for os=oshape_one:num_oshapes
    int_side_vbf_vals[os] = int_vs_side_vbf_vals_for_oshape(os, basis, vbf)
  end
  int_side_vbf_vals
end

function ref_side_vs_side_vbf_vals(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
  const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
  const side_side_vbf_vals = Array(Array{R,4}, num_oshapes) # by fe oshape, (side_face_1, monn_1, side_face_2, monn_2)
  for os=oshape_one:num_oshapes
    side_side_vbf_vals[os] = side_vs_side_vbf_vals_for_oshape(os, basis, vbf)
  end
  side_side_vbf_vals
end

# function reference_vbf_values(basis::WeakFunsPolyBasis, vbf::AbstractVariationalBilinearForm)
#   const int_int_vbf_vals = @spawn(TODO)
#   const vbf_symm = is_symmetric(vbf)
#   const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
#   const int_int_vbf_vals   = Array(Array{R,2}, num_oshapes) # vbf values indexed by fe oshape, (monn_1, monn_2)
#   const side_int_vbf_vals  = Array(Array{R,3}, num_oshapes) # by fe oshape, (side_face, side_monn, int_monn)
#   const int_side_vbf_vals  = Array(Array{R,3}, num_oshapes) # by fe oshape, (side_face, side_monn, int_monn) # [sic]
#   const side_side_vbf_vals = Array(Array{R,4}, num_oshapes) # by fe oshape, (side_face_1, monn_1, side_face_2, monn_2)
#   for os=oshape_one:num_oshapes
#     int_int_vbf_vals[os] = int_vs_int_vbf_vals_for_oshape(os, basis, vbf)
#     side_int_vbf_vals[os] = side_vs_int_vbf_vals_for_oshape(os, basis, vbf)
#     int_side_vbf_vals[os] = vbf_symm ? side_int_vbf_vals[os] : int_vs_side_vbf_vals_for_oshape(os, basis, vbf)
#     side_side_vbf_vals[os] = side_vs_side_vbf_vals_for_oshape(os, basis, vbf)
#   end
#   (int_int_vbf_vals, side_int_vbf_vals, int_side_vbf_vals, side_side_vbf_vals)
# end


# Returns an array of vbf values indexed by (monn_1, monn_2).
function int_vs_int_vbf_vals_for_oshape(fe_oshape::OShapeNum,
                                        basis::WeakFunsPolyBasis,
                                        vbf::AbstractVariationalBilinearForm)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const vbf_vals = Array(R, num_int_mons, num_int_mons)
  if is_symmetric(vbf)
    for monn_1=monnum(1):num_int_mons
      for monn_2=monnum(monn_1+1):num_int_mons
        const vbf_val = int_mon_vs_int_mon(fe_oshape, monn_1, monn_2, basis, vbf)
        vbf_vals[monn_1, monn_2] = vbf_val
        vbf_vals[monn_2, monn_1] = vbf_val
      end
      vbf_vals[monn_1, monn_1] = int_mon_vs_int_mon(fe_oshape, monn_1, monn_1, basis, vbf)
    end
  else
    for monn_1=monnum(1):num_int_mons, monn_2=monnum(1):num_int_mons
      vbf_vals[monn_1, monn_2] = int_mon_vs_int_mon(fe_oshape, monn_1, monn_2, basis, vbf)
    end
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face, side_monn, int_monn).
function side_vs_int_vbf_vals_for_oshape(fe_oshape::OShapeNum,
                                         basis::WeakFunsPolyBasis,
                                         vbf::AbstractVariationalBilinearForm)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)
  const vbf_vals = Array(R, num_sides, num_side_mons, num_int_mons)
  for int_monn=monnum(1):num_int_mons, sf=feface_one:num_sides, sf_monn=monnum(1):num_side_mons
    vbf_vals[sf, sf_monn, int_monn] = side_mon_vs_int_mon(fe_oshape, sf_monn, sf, int_monn, basis, vbf)
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face, side_monn, int_monn)
# The items are ordered so for substitutability with side_vs_int version above in the symmetric case.
function int_vs_side_vbf_vals_for_oshape(fe_oshape::OShapeNum,
                                         basis::WeakFunsPolyBasis,
                                         vbf::AbstractVariationalBilinearForm)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)
  const vbf_vals = Array(R, num_sides, num_side_mons, num_int_mons)
  for int_monn=monnum(1):num_int_mons, sf=feface_one:num_sides, sf_monn=monnum(1):num_side_mons
    vbf_vals[sf, sf_monn, int_monn] = int_mon_vs_side_mon(fe_oshape, int_monn, sf_monn, sf, basis, vbf)
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face_1, monn_1, side_face_2, monn_2).
function side_vs_side_vbf_vals_for_oshape(fe_oshape::OShapeNum,
                                          basis::WeakFunsPolyBasis,
                                          vbf::AbstractVariationalBilinearForm)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)

  const vbf_vals = Array(R, num_sides, num_side_mons, num_sides, num_side_mons)

  if !is_symmetric(vbf)
    for sf_1=feface_one:num_sides, sf_2=feface_one:num_sides,
        monn_1=monnum(1):num_side_mons, monn_2=monnum(1):num_side_mons
      vbf_vals[sf_1, monn_1, sf_2, monn_2] = side_mon_vs_side_mon(fe_oshape, monn_1, sf_1, monn_2, sf_2, basis, vbf)
    end
  else # vbf is symmetric
    # fill entries where sides differ (symmetric case)
    for sf_1=feface_one:num_sides, sf_2=fefacenum(sf_1+1):num_sides
      for monn_1=monnum(1):num_side_mons, monn_2=monnum(1):num_side_mons
        const vbf_val = side_mon_vs_side_mon(fe_oshape, monn_1, sf_1, monn_2, sf_2, basis, vbf)
        vbf_vals[sf_1, monn_1, sf_2, monn_2] = vbf_val
        vbf_vals[sf_2, monn_2, sf_1, monn_1] = vbf_val
      end
    end
    # fill entries where sides are the same (symmetric case)
    for sf=feface_one:num_sides
      for monn_1=monnum(1):num_side_mons
        for monn_2=monnum(monn_1+1):num_side_mons
          const vbf_val = side_mon_vs_side_mon(fe_oshape, monn_1, sf, monn_2, sf, basis, vbf)
          vbf_vals[sf, monn_1, sf, monn_2] = vbf_val
          vbf_vals[sf, monn_2, sf, monn_1] = vbf_val
        end
        vbf_vals[sf, monn_1, sf, monn_1] = side_mon_vs_side_mon(fe_oshape, monn_1, sf, monn_1, sf, basis, vbf)
      end
    end
  end

  vbf_vals
end

end # end of module
