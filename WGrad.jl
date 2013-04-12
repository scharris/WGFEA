module WGrad
export WGradSolver, wgrad

using Common
import Mesh, Mesh.AbstractMesh, Mesh.FERelFace
import Poly, Poly.Polynomial, Poly.Monomial, Poly.VectorMonomial, Poly.Nomial

# For a weak function v on a finite element T, the weak gradient of degree r
# of v on T is defined to be the polynomial vector function wgrad(v) in
# [P_r(T)]^d, such that:
#
# [WGRAD_DEF]
#   (wgrad(v), q)_T = -(v_0, div q)_T + <v_b, q.n>_bnd(T), for all q in [P_r(T)]^d
#
# By linearity of both sides in the right sides of the inner products, the above
# holds iff the same equation holds with the q's further restricted to be
# vector monomials on T (which are a monomial in one range component and 0 for
# all other components), which functions form a basis {b_i}_i for [P_r(T)]^d.
#
# Letting q = b_i, and writing wgrad(v) = sum_j eta_j b_j, we obtain a linear
# system from WGRAD_DEF which we can solve for the unkowns eta_j, with
# (b_i, b_j)_T being the matrix elements of the system, and the right hand side
# of WGRAD_DEF defining the right hand side column vector for the system.
#
# We will only actually need weak gradients of monomials or polynomials
# supported on a single face (interior or side) of a finite element. Thus v
# will be expressed below as a monomial or polynomial paired with the face of T
# on which it has the monomial or polynomial value.

type WGradSolver

  mesh::AbstractMesh

  test_vec_mons::Array{VectorMonomial,1}

  basis_vs_basis_ips::Matrix{R}

  test_vec_divs::Array{Polynomial,1}

  # Create a weak gradient solver for the specified weak gradient degree over the indicated mesh.
  function WGradSolver(wgrad_deg::Deg, mesh::AbstractMesh)
    const test_vmons = Poly.vector_monomials_of_degree_le(wgrad_deg, Mesh.space_dim(mesh))
    new(mesh,
        test_vmons,
        basis_vs_basis_ips(test_vmons, mesh),
        [Poly.divergence(q) for q in test_vmons])
  end
end

# Builds the basis vs. basis inner products matrix for the linear system described above.
function basis_vs_basis_ips(basis::Array{VectorMonomial,1}, mesh::AbstractMesh)
  const n = length(basis)
  const m = Array(R, n,n)
  for i=1:n
    const bi = basis[i]
    const bi_mon = bi.mon
    for j=i-1:-1:1
      const bj = basis[j]
      const ip = if Poly.monomials_in_same_component(bi, bj)
                   Mesh.integral_face_rel_x_face_rel_on_face(bi_mon, bj.mon, Mesh.interior_face, mesh)
                 else zeroR end
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = Mesh.integral_face_rel_x_face_rel_on_face(bi_mon, bi_mon, Mesh.interior_face, mesh)
  end
  m
end


# On the reference finite element for our mesh, obtains the weak gradient
# polynomial vector for the weak function v which is a monomial or polynomial
# on supporting face v_sface and 0 elsewhere.
function wgrad(v::Nomial, v_sface::FERelFace, solver::WGradSolver)
  const rhs = [wgrad_def_rhs_comp(v, v_sface, q_num, solver) for q_num in 1:length(solver.test_vec_mons)]
  # solve the linear system
  const sol_coefs = solver.basis_vs_basis_ips \ rhs
  Poly.linear_comb(sol_coefs, solver.test_vec_mons)
end

# Compute one component of the right hand side of the equation (WGRAD_DEF),
#           -(v_0, div q)_T + <v_b, q.n>_bnd(T),
# on the reference finite element of the mesh for weak function v. Here v
# is specified as a monomial or polynomial and the supporting face on which
# it takes the monomial or polynomial value.
function wgrad_def_rhs_comp(v::Nomial, v_sface::FERelFace, test_vec_mon_num::Integer, solver::WGradSolver)
  const q = solver.test_vec_mons[test_vec_mon_num]
  if v_sface == Mesh.interior_face
    # Interior supported v: only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
    const div_q = solver.test_vec_divs[test_vec_mon_num]
    -Mesh.integral_face_rel_x_face_rel_on_face(v, div_q, Mesh.interior_face, solver.mesh)
  else
    # Side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero.
    Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(v, q, v_sface, solver.mesh)
  end
end

end # end of module
