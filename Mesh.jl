module Mesh

export FENum, fe_num, no_fe,
       NBSideNum, nbside_num,
       FEFace, fe_face, interior_face, no_face,
       NBSideInclusions,
       AbstractMesh,
       space_dim,
       one_mon,
       num_fes,
       num_nb_sides,
       num_b_sides,
       num_side_faces_per_fe,
       fe_inclusions_of_nb_side,
       fe_inclusions_of_nb_side!,
       is_boundary_side,
       integral_on_ref_fe_face,
       integral_prod_on_fe_face,
       integral_prod_on_ref_fe_face,
       integral_on_ref_fe_side_vs_outward_normal,
       integral_prod_on_ref_fe_side_vs_outward_normal

using Common
import Poly.Monomial, Poly.Polynomial, Poly.VectorMonomial, Poly.PolynomialVector

# Number type for enumeration of all finite elements in a mesh.
typealias FENum Uint64
fe_num(i::Integer) = if i > 0 uint64(i) else error("finite element number out of range") end
const no_fe = zero(FENum)

# Number type for enumeration of all non-boundary sides in a mesh.
typealias NBSideNum Uint64
nbside_num(i::Integer) = if i > 0 uint64(i) else error("non-boundary side number out of range") end

# Number type for enumerating face positions (top, left, etc) on a finite element.
# For all meshes, the interior face is always numbered 0, while the sides are
# numbered 1 ... num_side_faces_per_fe(mesh).
typealias FEFace Uint8
fe_face(i::Integer) = if i >= 0 uint8(i) else error("face number out of range") end
const interior_face = zero(FEFace)
const no_face = 0xff

# Given a side not in the outside boundary, this structure represents the two finite elements
# which include the side, together with the side's face number in each of the finite elements.
type NBSideInclusions
  fe1::FENum
  face_in_fe1::FEFace
  fe2::FENum
  face_in_fe2::FEFace
end
NBSideInclusions() = NBSideInclusions(no_fe, no_face, no_fe, no_face)


# An abstract type representing any finite element mesh.
abstract AbstractMesh

# Generic functions which every specific mesh type should implement.

space_dim{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete)")

# The constantly-one monomial for the space dimension of the mesh.
one_mon{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete)")

num_fes{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete)")

num_nb_sides{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete")

num_b_sides{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete")

num_side_faces_per_fe{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete")

fe_inclusions_of_nb_side!{M <: AbstractMesh}(i::NBSideNum, mesh::M, nb_side_incls::NBSideInclusions) =
  error("not implemented, mesh implementation is incomplete")

is_boundary_side{M <: AbstractMesh}(fe::FENum, face::FEFace, mesh::M) =
  error("not implemented, mesh implementation is incomplete")

integral_on_ref_fe_face{M <: AbstractMesh}(m::Monomial, face::FEFace, mesh::M) =
  error("not implemented, mesh implementation is incomplete")

integral_prod_on_fe_face{M <: AbstractMesh}(f::Function, mon::Monomial, fe::FENum, face::FEFace, mesh::M) =
  error("not implemented, mesh implementation is incomplete")

integral_on_ref_fe_side_vs_outward_normal{M <: AbstractMesh}(vm::VectorMonomial, face::FEFace, mesh::M) =
  error("not implemented, mesh implementation is incomplete")



# The following have default implementations or are defined in terms of the required
# generic functions and usually won't need to be implemented for specific mesh types.

function integral_on_ref_fe_face{M <: AbstractMesh}(p::Polynomial, face::FEFace, mesh::M)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_on_ref_fe_face(p.mons[i], face, mesh)
  end
  sum
end

function integral_on_ref_fe_face{M <: AbstractMesh}(c::R, face::FEFace, mesh::M)
  if c == zeroR
    zeroR
  else
    c * integral_on_ref_fe_face(one_mon(mesh))
  end
end

# Implementations could specialize this to avoid creating the product monomial.
integral_prod_on_ref_fe_face{M <: AbstractMesh}(mon1::Monomial, mon2::Monomial, face::FEFace, mesh::M) =
  integral_on_ref_fe_face(mon1 * mon2, face, mesh)

function integral_prod_on_ref_fe_face{M <: AbstractMesh}(mon::Monomial, p::Polynomial, face::FEFace, mesh::M)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_prod_on_ref_fe_face(mon, p.mons[i], face, mesh)
  end
  sum
end

function integral_prod_on_ref_fe_face{M <: AbstractMesh}(p::Polynomial, mon::Monomial, face::FEFace, mesh::M)
  integral_prod_on_ref_fe_face(mon, p, face, mesh)
end

function integral_prod_on_ref_fe_face{M <: AbstractMesh}(p1::Polynomial, p2::Polynomial, face::FEFace, mesh::M)
  sum = zeroR
  for i=1:length(p1.coefs)
    sum += p1.coefs[i] * integral_prod_on_ref_fe_face(p1.mons[i], p2, face, mesh)
  end
  sum
end


# Functional (but memory allocating) variant of fe_inclusions_of_nb_side!.
fe_inclusions_of_nb_side{M <: AbstractMesh}(side_num::NBSideNum, mesh::M) =
  let side_incls = NBSideInclusions()
    fe_inclusions_of_nb_side!(side_num, mesh, side_incls)
    side_incls
  end

# For the next two functions, some subtype implementations (e.g. rectangular meshes) may be able to
# specialize for efficiency by returning 0 immediately in many cases without constructing v*q.

function integral_prod_on_ref_fe_side_vs_outward_normal(v::Monomial, q::VectorMonomial, side::FEFace, mesh::AbstractMesh)
  integral_on_ref_fe_side_vs_outward_normal(v * q, side, mesh)
end

function integral_prod_on_ref_fe_side_vs_outward_normal(v::Polynomial, q::VectorMonomial, side::FEFace, mesh::AbstractMesh)
  const prod_poly = q.mon * v
  const vm = VectorMonomial(prod_poly.mons[1], q.mon_pos)
  sum = prod_poly.coefs[1] * integral_on_ref_fe_side_vs_outward_normal(vm, side, mesh)
  for i=2:length(prod_poly.mons)
    vm.mon = prod_poly.mons[i]
    sum += prod_poly.coefs[i] * integral_on_ref_fe_side_vs_outward_normal(vm, side, mesh)
  end
  sum
end

end # end of module
