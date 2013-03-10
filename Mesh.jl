module Mesh

export FENum, fe_num, no_fe,
       SideNum, side_num,
       Face, no_face,
       NBSideInclusions,
       AbstractMesh,
       num_fes,
       num_nb_sides,
       fe_for_interior,
       fe_inclusions_of_nb_side,
       fe_inclusions_of_nb_side!,
       integral_on_ref_fe_interior,
       integral_on_ref_fe_side_vs_outward_normal

using Common
import Poly.Monomial, Poly.Polynomial, Poly.VectorMonomial

# Number type for enumeration of all finite elements in a mesh.
typealias FENum Uint64
fe_num(i::Integer) = uint64(i)
const no_fe = fe_num(0)

# Number type for enumeration of all (non-boundary) sides in a mesh.
typealias SideNum Uint64
side_num(i::Integer) = uint64(i)

# Number type for enumerating face positions (top, left, etc) on a single finite element.
typealias Face Uint8
const interior_face = 0x00
const no_face = 0xff

# Given a side not in the outside boundary, this structure represents the two finite elements
# which include the side, together with the side's face number in each of the finite elements.
type NBSideInclusions
  fe1::FENum
  face_in_fe1::Face
  fe2::FENum
  face_in_fe2::Face
end
NBSideInclusions() = NBSideInclusions(no_fe, no_face, no_fe, no_face)


# An abstract type representing any finite element mesh.
abstract AbstractMesh

# Generic functions which every specific mesh type should implement.

num_fes{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete)")
num_nb_sides{M <: AbstractMesh}(mesh::M) =
  error("not implemented, mesh implementation is incomplete")
fe_for_interior{M <: AbstractMesh}(intr_num::FENum, mesh::M) =
  error("not implemented, mesh implementation is incomplete")
fe_inclusions_of_nb_side!{M <: AbstractMesh}(i::SideNum, mesh::M, nb_side_incls::NBSideInclusions) =
  error("not implemented, mesh implementation is incomplete")
integral_on_ref_fe_interior{M <: AbstractMesh}(mon::Monomial, mesh::M) =
  error("not implemented, mesh implementation is incomplete")
integral_on_ref_fe_side_vs_outward_normal{M <: AbstractMesh}(vm::VectorMonomial, face::Face, mesh::M) =
  error("not implemented, mesh implementation is incomplete")

# The following are defined in terms of the required generic functions and usually won't need
# to be implemented for specific mesh types.

function integral_on_ref_fe_interior{M <: AbstractMesh}(p::Polynomial, mesh::M)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_on_ref_fe_interior(p.mons[i], mesh)
  end
  sum
end

# Functional (but memory allocating) variant of fe_inclusions_of_nb_side!.
fe_inclusions_of_nb_side{M <: AbstractMesh}(side_num::SideNum, mesh::M) =
  let side_incls = NBSideInclusions()
    fe_inclusions_of_nb_side!(side_num, mesh, side_incls)
    side_incls
  end


end # end of module
