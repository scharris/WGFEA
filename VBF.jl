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
       bel_vs_bel_transpose,
       bel_vs_bel_transpose_dense,
       make_interior_mon_side_projs

using Common
import Poly.Polynomial
import Mesh, Mesh.FENum, Mesh.FERelFace, Mesh.fe_face
import Proj
import WGBasis, WGBasis.BElNum, WGBasis.beln, WGBasis.WeakFunsPolyBasis, WGBasis.MonNum, WGBasis.mon_num

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
# bilinear form bf satisfying the Element Summability requirement.

function int_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  const fe1 = WGBasis.support_interior_num(bel_1, basis)
  const fe2 = WGBasis.support_interior_num(bel_2, basis)
  if fe1 != fe2
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let monn_1 = WGBasis.interior_mon_num(bel_1, basis)
        monn_2 = WGBasis.interior_mon_num(bel_2, basis)
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
  if side_face == Mesh.no_face # no fe includes both supports
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let side_monn = WGBasis.side_mon_num(sbel, basis)
        int_monn = WGBasis.interior_mon_num(ibel, basis)
      int_mon_vs_side_mon(int_num, int_monn, side_monn, side_face, basis, bf)
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
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let side_monn = WGBasis.side_mon_num(sbel, basis)
        int_monn = WGBasis.interior_mon_num(ibel, basis)
      side_mon_vs_int_mon(int_num, side_monn, side_face, int_monn, basis, bf)
    end
  end
end

function side_bel_vs_side_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, bf::BF)
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
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe1, basis, bf)
    elseif incls_1.fe1 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe1, monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe2, basis, bf)
    else
      zeroR
    end +
    # contribution from incls_1.fe2
    if incls_1.fe2 == incls_2.fe1
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe1, basis, bf)
    elseif incls_1.fe2 == incls_2.fe2
      side_mon_vs_side_mon(incls_1.fe2, monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe2, basis, bf)
    else
      zeroR
    end
  end
end

function bel_vs_bel_transpose{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, bf::BF)
  const num_int_bels = basis.num_interior_bels
  const first_nb_side_bel = basis.first_nb_side_bel

  # data arrays for construction of the sparse matrix
  const interacting_bel_pairs_ub = WGBasis.upper_bound_est_bel_pairs_supported_on_common_fe(basis)
  const row_nums = Array(Int, interacting_bel_pairs_ub)
  const col_nums = Array(Int, interacting_bel_pairs_ub)
  const nonzeros = Array(R, interacting_bel_pairs_ub)

  println(STDERR, "Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println(STDERR, "Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements, $(int64(basis.total_bels)^2) pairs.")
  println(STDERR, "Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")

  nnz = 0 # current number of non-zero values stored, the last stored position in the data arrays
  if !is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      if mod(nnz, 500) == 0 println(STDERR, "ii($i,$j)") end
      const vbf_ij = int_bel_vs_int_bel(beln(i), beln(j), basis, bf)
      if vbf_ij != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ij
        row_nums[nnz] = j
        col_nums[nnz] = i
      end
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      if mod(nnz, 500) == 0 println(STDERR, "si($i,$j)") end
      const bel_i = beln(i)
      const bel_j = beln(j)
      const vbf_ij = side_bel_vs_int_bel(bel_i, bel_j, basis, bf)
      if vbf_ij != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ij
        row_nums[nnz] = j
        col_nums[nnz] = i
      end
      const vbf_ji = int_bel_vs_side_bel(bel_j, bel_i, basis, bf)
      if vbf_ji != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ji
        row_nums[nnz] = i
        col_nums[nnz] = j
      end
    end
    for i=first_nb_side_bel:basis.total_bels, j=first_nb_side_bel:basis.total_bels
      if mod(nnz, 500) == 0 println(STDERR, "ss($i,$j)") end
      const vbf_ij = side_bel_vs_side_bel(beln(i), beln(j), basis, bf)
      if vbf_ij != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ij
        row_nums[nnz] = j
        col_nums[nnz] = i
      end
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_i = beln(i)
      for j=1:i-1
        if mod(nnz, 500) == 0 println(STDERR, "ii($i,$j)") end
        const vbf_ij = int_bel_vs_int_bel(bel_i, beln(j), basis, bf)
        if vbf_ij != zeroR
          nnz += 1
          nonzeros[nnz] = vbf_ij
          row_nums[nnz] = j
          col_nums[nnz] = i
          nnz += 1
          nonzeros[nnz] = vbf_ij
          row_nums[nnz] = i
          col_nums[nnz] = j
        end
      end
      const vbf_ii = int_bel_vs_int_bel(bel_i, bel_i, basis, bf)
      if vbf_ii != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ii
        row_nums[nnz] = i
        col_nums[nnz] = i
      end
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      if mod(nnz, 500) == 0 println(STDERR, "si($i,$j)") end
      const vbf_ij = side_bel_vs_int_bel(beln(i), beln(j), basis, bf)
      if vbf_ij != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ij
        row_nums[nnz] = j
        col_nums[nnz] = i
        nnz += 1
        nonzeros[nnz] = vbf_ij
        row_nums[nnz] = i
        col_nums[nnz] = j
      end
    end
    for i=first_nb_side_bel:basis.total_bels
      const bel_i = beln(i)
      for j=first_nb_side_bel:i-1
        if mod(nnz, 500) == 0 println(STDERR, "ss($i,$j)") end
        const vbf_ij = side_bel_vs_side_bel(bel_i, beln(j), basis, bf)
        if vbf_ij != zeroR
          nnz += 1
          nonzeros[nnz] = vbf_ij
          row_nums[nnz] = j
          col_nums[nnz] = i
          nnz += 1
          nonzeros[nnz] = vbf_ij
          row_nums[nnz] = i
          col_nums[nnz] = j
        end
      end
      const vbf_ii = side_bel_vs_side_bel(bel_i, bel_i, basis, bf)
      if vbf_ii != zeroR
        nnz += 1
        nonzeros[nnz] = vbf_ii
        row_nums[nnz] = i
        col_nums[nnz] = i
      end
    end
  end

  println(STDERR, "Done computing vbf bel vs bel matrix")

  sparse(row_nums[1:nnz],
         col_nums[1:nnz],
         nonzeros[1:nnz],
         basis.total_bels,
         basis.total_bels)
end


function bel_vs_bel_transpose_dense{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, bf::BF)
  const num_int_bels = basis.num_interior_bels
  const first_nb_side_bel = basis.first_nb_side_bel
  const m = Array(R, basis.total_bels, basis.total_bels)

  if !is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      const ip = int_bel_vs_int_bel(beln(i), beln(j), basis, bf)
      m[j,i] = ip
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const bel_i = beln(i)
      const bel_j = beln(j)
      m[j,i] = side_bel_vs_int_bel(bel_i, bel_j, basis, bf)
      m[i,j] = int_bel_vs_side_bel(bel_j, bel_i, basis, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=first_nb_side_bel:basis.total_bels
      m[j,i] = side_bel_vs_side_bel(beln(i), beln(j), basis, bf)
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_i = beln(i)
      for j=1:i-1
        const ip = int_bel_vs_int_bel(bel_i, beln(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = int_bel_vs_int_bel(bel_i, bel_i, basis, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const ip = side_bel_vs_int_bel(beln(i), beln(j), basis, bf)
      m[j,i] = ip
      m[i,j] = ip
    end
    for i=first_nb_side_bel:basis.total_bels
      const bel_i = beln(i)
      for j=first_nb_side_bel:i-1
        const ip = side_bel_vs_side_bel(bel_i, beln(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = side_bel_vs_side_bel(bel_i, bel_i, basis, bf)
    end
  end
  m
end


# TODO: unit tests
function make_interior_mon_side_projs(basis::WeakFunsPolyBasis)
  const mesh = basis.mesh
  const int_mons = WGBasis.interior_mons(basis)
  const num_int_mons = mon_num(length(int_mons))
  const projs_by_int_monn = Array(Array{Array{Polynomial,1},1}, num_int_mons)
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  for int_monn=mon_num(1):num_int_mons
    projs_by_oshape = Array(Array{Polynomial,1}, num_oshapes)
    for os=Mesh.oshape(1):num_oshapes
      const sides_per_fe = Mesh.num_side_faces_for_shape(os, mesh)
      const projs_by_side = Array(Polynomial, sides_per_fe)
      for sf=fe_face(1):sides_per_fe
        const side_mons = WGBasis.side_mons_for_oshape_side(os, sf, basis)
        const proj_coefs = Proj.project_interior_mon_onto_oshape_side(int_mons[int_monn], os, sf, basis)
        projs_by_side[sf] = Polynomial(side_mons, proj_coefs)
      end
      projs_by_oshape[os] = projs_by_side
    end
    projs_by_int_monn[int_monn] = projs_by_oshape
  end
  projs_by_int_monn
end

end # end of module
