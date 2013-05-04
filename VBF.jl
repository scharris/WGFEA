module VBF
export AbstractVariationalBilinearForm,
       is_symmetric,
       int_mon_vs_int_mon,
       side_mon_vs_int_mon,
       int_mon_vs_side_mon,
       side_mon_vs_side_mon,
       bel_vs_bel_transpose,
       bel_vs_bel_transpose_dense,
       make_interior_mon_side_projs

using Common
import Poly.Polynomial
import Mesh, Mesh.AbstractMesh, Mesh.OrientedShape, Mesh.FERelFace, Mesh.fe_face, Mesh.oshape, Mesh.fe_num
import Proj
import WGBasis, WGBasis.BElNum, WGBasis.beln, WGBasis.WeakFunsPolyBasis, WGBasis.MonNum, WGBasis.mon_num

abstract AbstractVariationalBilinearForm

# TODO: Rewrite this to document the even stronger assumption being used
#       in the code, which is that all the bf_T's are the *same* function
#       after translation, ie. they are all implementable as a single bf
#       on a reference element, so long as their input functions are translated
#       properly.
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

int_mon_vs_int_mon{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape,
                                                          monn_1::MonNum,
                                                          monn_2::MonNum,
                                                          basis::WeakFunsPolyBasis,
                                                          bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_int_mon{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           int_monn::MonNum,
                                                           basis::WeakFunsPolyBasis,
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

int_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape,
                                                           int_monn::MonNum,
                                                           side_monn::MonNum, side_face::FERelFace,
                                                           basis::WeakFunsPolyBasis,
                                                           bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

side_mon_vs_side_mon{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape,
                                                            monn_1::MonNum, side_face_1::FERelFace,
                                                            monn_2::MonNum, side_face_2::FERelFace,
                                                            basis::WeakFunsPolyBasis,
                                                            bf::BF) =
  error("not implemented, bilinear form implementation is incomplete")

#
##########################################################################


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


function bel_vs_bel_transpose{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, vbf::BF)
  # data arrays for construction of the sparse matrix
  const interacting_bel_pairs_ub = WGBasis.ub_estimate_num_bel_bel_common_support_fe_triplets(basis)
  const row_nums = Array(Int, interacting_bel_pairs_ub)
  const col_nums = Array(Int, interacting_bel_pairs_ub)
  const nonzeros = Array(R, interacting_bel_pairs_ub)

  println("Computing vbf bel vs bel matrix, with $(int64(basis.total_bels)) basis elements, $(int64(basis.total_bels)^2) pairs.")
  println("Mesh has $(int(basis.mesh.num_fes)) finite elements, and $(int(basis.mesh.num_nb_sides)) nb sides.")
  println("Basis has $(int(basis.mons_per_fe_interior)) monomials per interior, $(int(basis.mons_per_fe_side)) monomials per side.")
  println("Data arrays for non-zeros are of initial (upper bound) size $(int64(interacting_bel_pairs_ub)).")
  flush(STDOUT)

  const mesh = basis.mesh
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)

  # precompute vbf values for each fe oriented shape and pairing of monomial and support face with monomial and support face.
  const int_int_vbf_vals, side_int_vbf_vals, int_side_vbf_vals, side_side_vbf_vals = reference_vbf_values(basis, vbf)

  # work array for remembering which sides are nb sides within a finite element
  const is_nb_side = Array(Bool, max_num_shape_sides(mesh))

  nnz = 0 # current number of non-zero values stored, the last stored position in the data arrays

  for fe=fe_num(1):Mesh.num_fes(mesh)
    const fe_oshape = Mesh.oriented_shape_for_fe(fe, mesh)
    const fe_num_sides = Mesh.num_side_faces_for_shape(fe_oshape, mesh)

    # fill is_nb_side work array (up to this shapes # sides)
    for sf=fe_face(1):fe_num_sides  is_nb_side[sf] = !Mesh.is_boundary_side(fe, sf, mesh) end

    # fill interior vs interior matrix values
    for monn_1=mon_num(1):num_int_mons, monn_2=mon_num(1):num_int_mons
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
    for sf=fe_face(1):fe_num_sides if is_nb_side[sf]
      for sf_monn=mon_num(1):num_side_mons
        const sf_beln = WGBasis.side_mon_bel_num(fe, sf, sf_monn, basis)

        # fill side vs interior and interior vs side matrix values
        for int_monn=mon_num(1):num_int_mons
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
        for sf2=fe_face(1):fe_num_sides if is_nb_side[sf2]
          for sf2_monn=mon_num(1):num_side_mons
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

  println("Done computing vbf bel vs bel matrix")
  flush(STDOUT)

  sparse(row_nums[1:nnz],
         col_nums[1:nnz],
         nonzeros[1:nnz],
         basis.total_bels,
         basis.total_bels)
end


function reference_vbf_values{BF <: AbstractVariationalBilinearForm}(basis::WeakFunsPolyBasis, vbf::BF)
  const vbf_symm = is_symmetric(vbf)
  const num_oshapes = Mesh.num_oriented_element_shapes(basis.mesh)
  const int_int_vbf_vals   = Array(Array{R,2}, num_oshapes) # vbf values indexed by fe oshape, (monn_1, monn_2)
  const side_int_vbf_vals  = Array(Array{R,3}, num_oshapes) # vbf values indexed by fe oshape, (side_face, side_monn, int_monn)
  const int_side_vbf_vals  = Array(Array{R,3}, num_oshapes) # vbf values indexed by fe oshape, (side_face, side_monn, int_monn) # [sic]
  const side_side_vbf_vals = Array(Array{R,4}, num_oshapes) # vbf values indexed by fe oshape, (side_face_1, monn_1, side_face_2, monn_2)
  for os=oshape(1):num_oshapes
    int_int_vbf_vals[os] = int_vs_int_vbf_vals_for_oshape(os, basis, vbf)
    side_int_vbf_vals[os] = side_vs_int_vbf_vals_for_oshape(os, basis, vbf)
    int_side_vbf_vals[os] = vbf_symm ? side_int_vbf_vals[os] : int_vs_side_vbf_vals_for_oshape(os, basis, vbf)
    side_side_vbf_vals[os] = side_vs_side_vbf_vals_for_oshape(os, basis, vbf)
  end
  (int_int_vbf_vals, side_int_vbf_vals, int_side_vbf_vals, side_side_vbf_vals)
end


# Returns an array of vbf values indexed by (monn_1, monn_2).
function int_vs_int_vbf_vals_for_oshape{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape, basis::WeakFunsPolyBasis, vbf::BF)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const vbf_vals = Array(R, num_int_mons, num_int_mons)
  if is_symmetric(vbf)
    for monn_1=mon_num(1):num_int_mons
      for monn_2=mon_num(monn_1+1):num_int_mons
        const vbf_val = int_mon_vs_int_mon(fe_oshape, monn_1, monn_2, basis, vbf)
        vbf_vals[monn_1, monn_2] = vbf_val
        vbf_vals[monn_2, monn_1] = vbf_val
      end
      vbf_vals[monn_1, monn_1] = int_mon_vs_int_mon(fe_oshape, monn_1, monn_1, basis, vbf)
    end
  else
    for monn_1=mon_num(1):num_int_mons, monn_2=mon_num(1):num_int_mons
      vbf_vals[monn_1, monn_2] = int_mon_vs_int_mon(fe_oshape, monn_1, monn_2, basis, vbf)
    end
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face, side_monn, int_monn).
function side_vs_int_vbf_vals_for_oshape{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape, basis::WeakFunsPolyBasis, vbf::BF)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)
  const vbf_vals = Array(R, num_sides, num_side_mons, num_int_mons)
  for int_monn=mon_num(1):num_int_mons, sf=fe_face(1):num_sides, sf_monn=mon_num(1):num_side_mons
    vbf_vals[sf, sf_monn, int_monn] = side_mon_vs_int_mon(fe_oshape, sf_monn, sf, int_monn, basis, vbf)
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face, side_monn, int_monn) [ordered so for substitutability with side_vs_int version above].
function int_vs_side_vbf_vals_for_oshape{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape, basis::WeakFunsPolyBasis, vbf::BF)
  const num_int_mons = WGBasis.mons_per_fe_interior(basis)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)
  const vbf_vals = Array(R, num_sides, num_side_mons, num_int_mons)
  for int_monn=mon_num(1):num_int_mons, sf=fe_face(1):num_sides, sf_monn=mon_num(1):num_side_mons
    vbf_vals[sf, sf_monn, int_monn] = int_mon_vs_side_mon(fe_oshape, int_monn, sf_monn, sf, basis, vbf)
  end
  vbf_vals
end

# Returns an array of vbf values indexed by (side_face_1, monn_1, side_face_2, monn_2).
function side_vs_side_vbf_vals_for_oshape{BF <: AbstractVariationalBilinearForm}(fe_oshape::OrientedShape, basis::WeakFunsPolyBasis, vbf::BF)
  const num_side_mons = WGBasis.mons_per_fe_side(basis)
  const num_sides = Mesh.num_side_faces_for_shape(fe_oshape, basis.mesh)

  const vbf_vals = Array(R, num_sides, num_side_mons, num_sides, num_side_mons)

  if !is_symmetric(vbf)
    for sf_1=fe_face(1):num_sides, sf_2=fe_face(1):num_sides,
        monn_1=mon_num(1):num_side_mons, monn_2=mon_num(1):num_side_mons
      vbf_vals[sf_1, monn_1, sf_2, monn_2] = side_mon_vs_side_mon(fe_oshape, monn_1, sf_1, monn_2, sf_2, basis, vbf)
    end
  else # vbf is symmetric
    # fill entries where sides differ (symmetric case)
    for sf_1=fe_face(1):num_sides, sf_2=fe_face(sf_1+1):num_sides
      for monn_1=mon_num(1):num_side_mons, monn_2=mon_num(1):num_side_mons
        const vbf_val = side_mon_vs_side_mon(fe_oshape, monn_1, sf_1, monn_2, sf_2, basis, vbf)
        vbf_vals[sf_1, monn_1, sf_2, monn_2] = vbf_val
        vbf_vals[sf_2, monn_2, sf_1, monn_1] = vbf_val
      end
    end
    # fill entries where sides are the same (symmetric case)
    for sf=fe_face(1):num_sides
      for monn_1=mon_num(1):num_side_mons
        for monn_2=mon_num(monn_1+1):num_side_mons
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


function max_num_shape_sides(mesh::AbstractMesh)
  max_sides= 0
  for os=oshape(1):Mesh.num_oriented_element_shapes(mesh)
    max_sides = max(max_sides, Mesh.num_side_faces_for_shape(os, mesh))
  end
  max_sides
end



# The remaining functions are used only to implement a dense matrix implementation of the vbf bel vs bel matrix which uses
# a global bel vs. bel comparison, as opposed to the more efficient fe-wise computations of the new sparse approach.  It is
# usefull mainly to test the fe-wise sparse implementation.  These functions could usefully be moved into a test file.


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


function int_bel_vs_int_bel{BF <: AbstractVariationalBilinearForm}(bel_1::BElNum, bel_2::BElNum, basis::WeakFunsPolyBasis, bf::BF)
  const fe1 = WGBasis.support_interior_num(bel_1, basis)
  const fe2 = WGBasis.support_interior_num(bel_2, basis)
  if fe1 != fe2
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let monn_1 = WGBasis.interior_mon_num(bel_1, basis)
        monn_2 = WGBasis.interior_mon_num(bel_2, basis)
      int_mon_vs_int_mon(Mesh.oriented_shape_for_fe(fe1,basis.mesh), monn_1, monn_2, basis, bf)
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
      int_mon_vs_side_mon(Mesh.oriented_shape_for_fe(int_num,basis.mesh), int_monn, side_monn, side_face, basis, bf)
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
      side_mon_vs_int_mon(Mesh.oriented_shape_for_fe(int_num,basis.mesh), side_monn, side_face, int_monn, basis, bf)
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
      side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe1,basis.mesh), monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe1, basis, bf)
    elseif incls_1.fe1 == incls_2.fe2
      side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe1,basis.mesh), monn_1, incls_1.face_in_fe1, monn_2, incls_2.face_in_fe2, basis, bf)
    else
      zeroR
    end +
    # contribution from incls_1.fe2
    if incls_1.fe2 == incls_2.fe1
      side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe2,basis.mesh), monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe1, basis, bf)
    elseif incls_1.fe2 == incls_2.fe2
      side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe2,basis.mesh), monn_1, incls_1.face_in_fe2, monn_2, incls_2.face_in_fe2, basis, bf)
    else
      zeroR
    end
  end
end


end # end of module
