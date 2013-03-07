module WGrad
export WGradSolver, wgrad

using RMesh
using Poly

# For a weak function b on a finite element T, the weak gradient (of degree r)
# of b on T is defined to be the polynomial vector function wgrad b in [P_r(T)]^d,
# such that:
# [WGRAD_DEF]
#   (wgrad b, q)_T = -(b_0, div q)_T + <b_b, q.n>_bnd(T),
# for all polynomial vectors q in [P_r(T)]^d.


type WGradSolver

  mesh::RectMesh

  test_vec_mons::Array{VectorMonomial,1}

  basis_vs_basis_ips::Matrix{R}

  test_vec_divs::Array{Polynomial,1}

  # Create a weak gradient solver for the specified weak gradient degree and domain dimension.
  # In the WG method, wgrad_deg = k-1 where k is the max degree of the approximating polynomials
  # on finite element interiors.
  function WGradSolver(wgrad_deg::Pwr, dom_dim::Dim, mesh::RectMesh)
    vmons = vector_monomials_of_degree_le(wgrad_deg, dom_dim)
    basis_ips = basis_vs_basis_ips(vmons, mesh)
    divs = [divergence(q) for q in vmons]

    new(mesh,
        vmons,
        basis_ips,
        divs)
  end
end

function basis_vs_basis_ips(basis::Array{VectorMonomial,1}, mesh::RectMesh)
  const n = length(basis)
  const m = Array(R, n,n)
  for i=1:n
    for j=i-1:-1:1
      const ip = integrate_on_llo_rect(dot(basis[i], basis[j]), mesh.fe_width, mesh.fe_height)
      m[i,j] = ip
      m[j,i] = ip
    end
    m[i,i] = integrate_on_llo_rect(dot(basis[i], basis[i]), mesh.fe_width, mesh.fe_height)
  end
  m
end

function wgrad(b_mon::Monomial, b_sface::FaceNum, wgrad_solver::WGradSolver)
  const fe_w = wgrad_solver.mesh.fe_width
  const fe_h = wgrad_solver.mesh.fe_height
  const rhs = [ip_wgrad_bel_with_test_vec_mon(b_mon, b_sface, q_num, wgrad_solver, fe_w, fe_h)
                for q_num in 1:length(wgrad_solver.test_vec_mons)]
  # solve the linear system
  const coefs = wgrad_solver.basis_vs_basis_ips \ rhs
  # apply coeficients to the vector monomials basis functions
  linear_comb(coefs, wgrad_solver.test_vec_mons)
end

# Computes the right hand side of the equation (WGRAD_DEF).
# Here the weak function is a main basis element b and is described as a monomial paired
# with the supporting face on which it takes this monomial definition (being 0 elsewhere).
# Both monomials take their domain values relative to the lower left corner of the finite
# element (which is why only the finite element's width and height are needed).
function ip_wgrad_bel_with_test_vec_mon(b_mon::Monomial, b_sface::FaceNum,
                                        vec_mon_num::Int,
                                        wgrad_solver::WGradSolver,
                                        fe_w::R, fe_h::R)
  const q = wgrad_solver.test_vec_mons[vec_mon_num]

  # For interior supported b, only the -(b_0, div q)_T term is non-zero in the rhs of (WGRAD_DEF).
  if b_sface == RMesh.interior_face
    const div_q = wgrad_solver.test_vec_divs[vec_mon_num]
    -integrate_on_llo_rect(b_mon * div_q, fe_w, fe_h)

  # only the <b_b, q.n>_bnd(T) term is non-zero for the remaining (edge supported b) cases

  elseif b_sface == RMesh.top_face
    # y = fe_h on the top face, and q.n = q[2] there.
    dim_red_intgd = reduce_dim_by_fixing(dim(2), fe_h, b_mon * q[2])
    integrate_on_llo_rect(dim_red_intgd, fe_w) # integrate over the remaining non-fixed dimensions

  elseif b_sface == RMesh.right_face
    # x=fe_w on the right face, and q.n = q[1] there.
    dim_red_intgd = reduce_dim_by_fixing(dim(1), fe_w, b_mon * q[1])
    integrate_on_llo_rect(dim_red_intgd, fe_h) # integrate over the remaining non-fixed dimensions

  elseif b_sface == RMesh.bottom_face
    # y=0 on the bottom face, and q.n = -q[2] there.
    # optimization (to avoid allocations): if there is a non-zero exponent for y in either monomial factor then the integral is 0.
    if b_mon.exps[2] != 0 || let qcomp = q[2]; length(qcomp.mons) == 1 && qcomp.mons[1].exps[2] != 0 end
      zeroR
    else
      dim_red_intgd = reduce_dim_by_fixing(dim(2), zeroR, -(b_mon * q[2]))
      integrate_on_llo_rect(dim_red_intgd, fe_w) # integrate over the remaining non-fixed dimensions
    end

  elseif b_sface == RMesh.left_face
    # x=0 on the left face, and q.n = -q[1] there.
    # optimization (to avoid allocations): if there is a non-zero exponent for x in either monomial factor then the integral is 0.
    if b_mon.exps[1] != 0 || let qcomp = q[1]; length(qcomp.mons) == 1 && qcomp.mons[1].exps[1] != 0 end
      zeroR
    else
      dim_red_intgd = reduce_dim_by_fixing(dim(1), zeroR, -(b_mon * q[1]))
      integrate_on_llo_rect(dim_red_intgd, fe_h) # integrate over the remaining non-fixed dimensions
    end

  else
    error("invalid face")
  end
end

end # end of module
