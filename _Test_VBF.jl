using Base.Test
using Common
import Mesh
import RMesh.RectMesh, RMesh.mesh_ldims
import WGBasis, WGBasis.BElNum, WGBasis.belnum, WGBasis.WeakFunsPolyBasis
import VBF, VBF.AbstractVariationalBilinearForm
import VBF_a_s.a_s


# The following function and supporting functions are used for testing the sparse matrix creation in the VBF module,
# which is to be tested against the simpler dense implementation here.

function bel_vs_bel_transpose_dense(basis::WeakFunsPolyBasis, bf::AbstractVariationalBilinearForm)
  const num_int_bels = basis.num_interior_bels
  const first_nb_side_bel = basis.first_nb_side_bel
  const m = Array(R, basis.total_bels, basis.total_bels)

  if !VBF.is_symmetric(bf)
    for i=1:num_int_bels, j=1:num_int_bels
      const ip = int_bel_vs_int_bel(belnum(i), belnum(j), basis, bf)
      m[j,i] = ip
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const bel_i = belnum(i)
      const bel_j = belnum(j)
      m[j,i] = side_bel_vs_int_bel(bel_i, bel_j, basis, bf)
      m[i,j] = int_bel_vs_side_bel(bel_j, bel_i, basis, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=first_nb_side_bel:basis.total_bels
      m[j,i] = side_bel_vs_side_bel(belnum(i), belnum(j), basis, bf)
    end
  else # bf is symmetric
    for i=1:num_int_bels
      const bel_i = belnum(i)
      for j=1:i-1
        const ip = int_bel_vs_int_bel(bel_i, belnum(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = int_bel_vs_int_bel(bel_i, bel_i, basis, bf)
    end
    for i=first_nb_side_bel:basis.total_bels, j=1:num_int_bels
      const ip = side_bel_vs_int_bel(belnum(i), belnum(j), basis, bf)
      m[j,i] = ip
      m[i,j] = ip
    end
    for i=first_nb_side_bel:basis.total_bels
      const bel_i = belnum(i)
      for j=first_nb_side_bel:i-1
        const ip = side_bel_vs_side_bel(bel_i, belnum(j), basis, bf)
        m[j,i] = ip
        m[i,j] = ip
      end
      m[i,i] = side_bel_vs_side_bel(bel_i, bel_i, basis, bf)
    end
  end
  m
end


function int_bel_vs_int_bel(bel_1::BElNum,
                            bel_2::BElNum,
                            basis::WeakFunsPolyBasis,
                            bf::AbstractVariationalBilinearForm)
  const fe1 = WGBasis.support_interior_num(bel_1, basis)
  const fe2 = WGBasis.support_interior_num(bel_2, basis)
  if fe1 != fe2
    zeroR # by locality property
  else
    # By common support summability property, we only need the contribution from the single common fe.
    let monn_1 = WGBasis.interior_mon_num(bel_1, basis)
        monn_2 = WGBasis.interior_mon_num(bel_2, basis)
      VBF.int_mon_vs_int_mon(Mesh.oriented_shape_for_fe(fe1,basis.mesh), monn_1, monn_2, basis, bf)
    end
  end
end


function int_bel_vs_side_bel(ibel::BElNum,
                             sbel::BElNum,
                             basis::WeakFunsPolyBasis,
                             bf::AbstractVariationalBilinearForm)
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
      VBF.int_mon_vs_side_mon(Mesh.oriented_shape_for_fe(int_num,basis.mesh), int_monn, side_monn, side_face, basis, bf)
    end
  end
end


function side_bel_vs_int_bel(sbel::BElNum,
                             ibel::BElNum,
                             basis::WeakFunsPolyBasis,
                             bf::AbstractVariationalBilinearForm)
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
      VBF.side_mon_vs_int_mon(Mesh.oriented_shape_for_fe(int_num,basis.mesh), side_monn, side_face, int_monn, basis, bf)
    end
  end
end


function side_bel_vs_side_bel(bel_1::BElNum,
                              bel_2::BElNum,
                              basis::WeakFunsPolyBasis,
                              bf::AbstractVariationalBilinearForm)
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
      VBF.side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe1,basis.mesh),
                               monn_1, incls_1.face_in_fe1,
                               monn_2, incls_2.face_in_fe1,
                               basis,
                               bf)
    elseif incls_1.fe1 == incls_2.fe2
      VBF.side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe1,basis.mesh),
                               monn_1, incls_1.face_in_fe1,
                               monn_2, incls_2.face_in_fe2,
                               basis,
                               bf)
    else
      zeroR
    end +
    # contribution from incls_1.fe2
    if incls_1.fe2 == incls_2.fe1
      VBF.side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe2,basis.mesh),
                                                          monn_1, incls_1.face_in_fe2,
                                                          monn_2, incls_2.face_in_fe1,
                                                          basis,
                                                          bf)
    elseif incls_1.fe2 == incls_2.fe2
      VBF.side_mon_vs_side_mon(Mesh.oriented_shape_for_fe(incls_1.fe2,basis.mesh),
                                                          monn_1, incls_1.face_in_fe2,
                                                          monn_2, incls_2.face_in_fe2,
                                                          basis,
                                                          bf)
    else
      zeroR
    end
  end
end


mesh = RectMesh([0.,0.], [1.,1.], mesh_ldims(5,5))
basis = WeakFunsPolyBasis(deg(3), deg(2), mesh)
vbf = a_s(basis)
sparse_m = VBF.bel_vs_bel_transpose(basis, vbf)
dense_m = bel_vs_bel_transpose_dense(basis, vbf)

# println("Max abs diff: ", max(abs(sparse_m - dense_m)))
@test max(abs(sparse_m - dense_m)) < 1e-12
