module WGrad
export WGradSolver, wgrad

using Common
import Mesh, Mesh.AbstractMesh, Mesh.FEFace
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
# polynomial vector for the weak function v which is a monomial or polynomial
# on supporting face v_sface and 0 elsewhere.
function wgrad(v::Nomial, v_sface::FEFace, solver::WGradSolver)
  const rhs = [wgrad_def_rhs(v, v_sface, q_num, solver) for q_num in 1:length(solver.test_vec_mons)]
  # solve the linear system
  const sol_coefs = solver.basis_vs_basis_ips \ rhs
  Poly.linear_comb(sol_coefs, solver.test_vec_mons)
end

# Computes the right hand side of the equation (WGRAD_DEF) on the reference
# finite element of the mesh for a weak function v. Here v is specified as
# a monomial or polynomial and the supporting face on which it takes the
# monomial or polynomial value.
function wgrad_def_rhs(v::Nomial, v_sface::FEFace,
                       test_vec_mon_num::Int,
                       solver::WGradSolver)
  const q = solver.test_vec_mons[test_vec_mon_num]

  # For interior supported v, only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
  if v_sface == Mesh.interior_face
    const div_q = solver.test_vec_divs[test_vec_mon_num]
    -Mesh.integral_on_ref_fe_interior(v * div_q, solver.mesh)
  else
    # For side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero.
    Mesh.integral_prod_on_ref_fe_side_vs_outward_normal(v, q, v_sface, solver.mesh)
  end
end

end # end of module
