module WGrad
export WGradSolver, wgrad

using Common
import Mesh, Mesh.AbstractMesh, Mesh.RelFace, Mesh.OrientedShape
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

immutable WGradSolver

  mesh::AbstractMesh

  basis_vmons::Array{VectorMonomial,1}

  ips_basis_vs_basis::Array{Matrix{R},1} # indexed by fe oriented shape

  basis_divs::Array{Polynomial,1}

  # Create a weak gradient solver for the specified weak gradient degree over the indicated mesh.
  function WGradSolver(wgrad_deg::Deg, mesh::AbstractMesh)
    const basis_vmons = Poly.vector_mons_of_deg_le(wgrad_deg, Mesh.space_dim(mesh))
    new(mesh,
        basis_vmons,
        make_vmon_ips_by_oshape(basis_vmons, mesh),
        [Poly.divergence(q) for q in basis_vmons])
  end
end

# Builds the array of vmon vs. vmon inner products matrices for the linear system described above, indexed by fe oriented shape.
function make_vmon_ips_by_oshape(vmons::Array{VectorMonomial,1}, mesh::AbstractMesh)
  const num_oshapes = Mesh.num_oriented_element_shapes(mesh)
  const ips_by_oshape = Array(Matrix{R}, num_oshapes)
  const num_vmons = length(vmons)
  for os=Mesh.oshape(1):num_oshapes
    const m = Array(R, num_vmons,num_vmons)
    for i=1:num_vmons
      const bi = vmons[i]
      const bi_mon = bi.mon
      for j=i-1:-1:1
        const bj = vmons[j]
        const ip = if Poly.mons_in_same_comp_of_vmons(bi, bj)
                     Mesh.integral_face_rel_x_face_rel_on_face(bi_mon, bj.mon, os, Mesh.interior_face, mesh)
                   else zeroR end
        m[i,j] = ip
        m[j,i] = ip
      end
      m[i,i] = Mesh.integral_face_rel_x_face_rel_on_face(bi_mon, bi_mon, os, Mesh.interior_face, mesh)
    end
    ips_by_oshape[os] = m
  end
  ips_by_oshape
end


# On the reference finite element for our mesh, obtains the weak gradient
# polynomial vector for the weak function v which is a monomial or polynomial
# on supporting face v_sface and 0 elsewhere.
function wgrad(v::Nomial, fe_oshape::OrientedShape, v_sface::RelFace, solver::WGradSolver)
  const rhs = [wgrad_def_rhs_comp(v, fe_oshape, v_sface, vmon_num, solver)
               for vmon_num in 1:length(solver.basis_vmons)]
  # solve the linear system
  const sol_coefs = solver.ips_basis_vs_basis[fe_oshape] \ rhs
  Poly.linear_comb(sol_coefs, solver.basis_vmons)
end

# Compute one component of the right hand side of the equation (WGRAD_DEF),
#           -(v_0, div q)_T + <v_b, q.n>_bnd(T),
# on the reference finite element of the mesh for weak function v. Here v
# is specified as a monomial or polynomial and the supporting face on which
# it takes the monomial or polynomial value.
function wgrad_def_rhs_comp(v::Nomial, fe_oshape::OrientedShape, v_sface::RelFace, vmon_num::Integer, solver::WGradSolver)
  const q = solver.basis_vmons[vmon_num]
  if v_sface == Mesh.interior_face
    # Interior supported v: only the -(v_0, div q)_T term can be non-zero in the rhs of (WGRAD_DEF).
    const div_q = solver.basis_divs[vmon_num]
    -Mesh.integral_face_rel_x_face_rel_on_face(v, div_q, fe_oshape, Mesh.interior_face, solver.mesh)
  else
    # Side supported v: only the <v_b, q.n>_bnd(T) term can be non-zero.
    Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_side(v, q, fe_oshape, v_sface, solver.mesh)
  end
end

end # end of module
