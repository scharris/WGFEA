module WGrad
export WGradSolver, wgrad

using Common
import Mesh, Mesh.AbstractMesh, Mesh.Face
import Poly, Poly.Polynomial, Poly.Monomial, Poly.VectorMonomial

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
# We will only actually need weak gradients of basis functions for the
# appoximation space V_h, which basis functions will be a monomial on one face
# (interior or side) of T and 0 on other faces. Thus v will be expressed below
# as a monomial paired with the face of T on which it has the monomial value.

type WGradSolver

  mesh::AbstractMesh

  test_vec_mons::Array{VectorMonomial,1}

  basis_vs_basis_ips::Matrix{R}

  test_vec_divs::Array{Polynomial,1}

  # Create a weak gradient solver for the specified weak gradient degree and
  # domain dimension. In the WG method, wgrad_deg = k-1 where k is the max
  # degree of the approximating polynomials on finite element interiors.
  function WGradSolver(wgrad_deg::Deg, dom_dim::Dim, mesh::AbstractMesh)
    const vms = Poly.vector_monomials_of_degree_le(wgrad_deg, dom_dim)
    new(mesh,
        vms,
        basis_vs_basis_ips(vms, mesh),
        [Poly.divergence(q) for q in vms])
  end
end

# Builds the matrix for the linear system described above.
function basis_vs_basis_ips(basis::Array{VectorMonomial,1}, mesh::AbstractMesh)
  const n = length(basis)
  const m = Array(R, n,n)
  for i=1:n
    const bi = basis[i]
    for j=i-1:-1:1
      const bj = basis[j]
      const ip = if Poly.monomials_in_same_component(bi, bj)
                   Mesh.integral_on_ref_fe_interior(bi.mon * bj.mon, mesh)
                 else zeroR end
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = Mesh.integral_on_ref_fe_interior(bi.mon * bi.mon, mesh)
  end
  m
end

# On the reference finite element for our mesh, obtains the weak gradient
# polynomial vector for the weak function which is a monomial v_mon on
# supporting face v_sface and 0 elsewhere.
function wgrad(v_mon::Monomial, v_sface::Face, wgrad_solver::WGradSolver)
  const rhs = [ip_wgrad_bel_with_test_vec_mon(v_mon, v_sface, q_num, wgrad_solver)
                for q_num in 1:length(wgrad_solver.test_vec_mons)]
  # solve the linear system
  const coefs = wgrad_solver.basis_vs_basis_ips \ rhs

  Poly.linear_comb(coefs, wgrad_solver.test_vec_mons)
end

# Computes the right hand side of the equation (WGRAD_DEF) on the reference
# finite element of the mesh. Here the weak function is a main basis element v
# and is described as a monomial paired with the supporting face on which it
# takes this monomial definition, being 0 elsewhere.
function ip_wgrad_bel_with_test_vec_mon(v_mon::Monomial, v_sface::Face,
                                        test_vec_mon_num::Int,
                                        wgrad_solver::WGradSolver)
  const q = wgrad_solver.test_vec_mons[test_vec_mon_num]

  # For interior supported v, only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
  if v_sface == Mesh.interior_face
    const div_q = wgrad_solver.test_vec_divs[test_vec_mon_num]
    -Mesh.integral_on_ref_fe_interior(v_mon * div_q, wgrad_solver.mesh)
  else
    # For side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero.
    Mesh.integral_on_ref_fe_side_vs_outward_normal(v_mon * q, v_sface, wgrad_solver.mesh)
  end
end

end # end of module
