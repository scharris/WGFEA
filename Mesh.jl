module Mesh

export FENum, fe_num, no_fe,
       NBSideNum, nb_side_num,
       FERelFace, fe_face, interior_face, no_face,
       OrientedShape, oshape,
       NBSideInclusions,
       AbstractMesh,
       space_dim,
       one_mon,
       num_fes,
       fe_interior_origin,
       fe_interior_origin!,
       shape_diameter_inv,
       max_fe_diameter,
       num_nb_sides,
       num_oriented_element_shapes,
       oriented_shape_for_fe,
       num_side_faces_for_fe,
       num_side_faces_for_shape,
       max_num_shape_sides,
       dependent_dim_for_nb_side,
       dependent_dim_for_oshape_side,
       fe_inclusions_of_nb_side,
       nb_side_num_for_fe_side,
       is_boundary_side,
       num_boundary_sides,
       num_non_boundary_sides_for_fe,
       integral_face_rel_on_face,
       integral_global_x_face_rel_on_fe_face,
       integral_face_rel_x_face_rel_on_face,
       integral_side_rel_x_fe_rel_vs_outward_normal_on_side,
       integral_fe_rel_x_side_rel_on_side


using Common
import Poly.Monomial, Poly.Polynomial, Poly.VectorMonomial, Poly.PolynomialVector, Poly.NomialOrConst

# Number type for enumeration of all finite elements in a mesh.
typealias FENum Uint64
fe_num(i::Integer) = if i > 0 convert(Uint64, i) else error("finite element number out of range") end
const no_fe = zero(FENum)

# Number type for enumeration of all non-boundary sides in a mesh.
typealias NBSideNum Uint64
no_nb_side = uint64(0)
nb_side_num(i::Integer) = if i > 0 convert(Uint64, i) else error("non-boundary side number out of range") end

# Finite Element Relative Face - Number type for enumerating the interior and
# side parts of a finite element relative to the element itself. For all
# meshes, the interior face is always numbered 0, while the sides are numbered
# 1 ... num_side_faces_for_fe(fe, mesh) for a given fe.
typealias FERelFace Uint8
fe_face(i::Integer) = if i >= 0 convert(Uint8, i) else error("face number out of range") end
const interior_face = zero(FERelFace)
const no_face = 0xff

typealias OrientedShape Int8
oshape(i::Integer) = if i > 127 || i <= 0 error("oriented shape number out of range") else convert(Int8, i) end

# Given a side not in the outside boundary, this structure represents the two
# finite elements which include the side, together with the side's face number
# in each of the finite elements.
immutable NBSideInclusions
  nb_side_num::NBSideNum
  fe1::FENum
  face_in_fe1::FERelFace
  fe2::FENum
  face_in_fe2::FERelFace
end
NBSideInclusions() = NBSideInclusions(no_nb_side, no_fe, no_face, no_fe, no_face)


# An abstract type representing any finite element mesh.
abstract AbstractMesh

####################################################################
# Generic functions which every specific mesh type should implement.

space_dim(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

# The constantly-one monomial for the space dimension of the mesh.
one_mon(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

num_fes(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

num_nb_sides(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

num_oriented_element_shapes(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

oriented_shape_for_fe(fe::FENum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

#num_fes_of_oriented_shape(oshape::OrientedShape, mesh::AbstractMesh) =
#  error("not implemented, mesh implementation is incomplete")

num_side_faces_for_fe(fe::FENum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

num_side_faces_for_shape(oshape::OrientedShape, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")


# The dependent_dim* functions provide, for a given side, the dimension j which is
# affine-dependent on the other dimensions on the side. That is, the function
# returns a j for which c_0,...,c_d exist such that
#             x_j = c_0 + sum_{i=1..d} c_i x_i for all (x_1,...,x_d) in the side.
# There may be more than one such coordinate number, in which case any one of
# these is returned.
dependent_dim_for_nb_side(i::NBSideNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

dependent_dim_for_oshape_side(fe_oshape::OrientedShape, side_face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

fe_inclusions_of_nb_side(side_num::NBSideNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Return non-boundary side number of the indicated fe relative side, or 0 if the side is a boundary side.
nb_side_num_for_fe_side(fe::FENum, side_face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")


is_boundary_side(fe::FENum, face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

num_boundary_sides(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

shape_diameter_inv(oshape::OrientedShape, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

max_fe_diameter(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

fe_interior_origin!(fe::FENum, fill::Vector{R}, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete)")

# Integration Functions

# Integrals are named according to the way in which their function arguments
# interpret their input values. Face relative arguments are interpreted with
# the origin for their input values being chosen by the mesh as a function of
# the absolute face, meaning the actual side or interior and not just its
# representation as a relative face in a finite element. That is, face relative
# arguments are implicitly pre-composed with a translation x -> x - o(F), where
# o(F) is a point which is chosen by the mesh for the face F, usually most
# naturally within the face itself. Likewise finite element relative functions
# employ an origin which is chosen for the finite element as a whole by the
# mesh.  Global functions are not pre-composed with any translation function.

# Integrate a monomial on the indicated face of the reference element, with the
# monomial interpreted locally on the face.
integral_face_rel_on_face(m::Monomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a global function f on the indicated finite element face, multiplied against
# a monomial interpreted locally on the face.
integral_global_x_face_rel_on_fe_face(f::Function, mon::Monomial, fe::FENum, face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a side-local monomial vs. a vector monomial interpreted relative to the
# entire finite element, dot multiplied with the outward normal for the side.
integral_side_rel_x_fe_rel_vs_outward_normal_on_side(v::Monomial, q::VectorMonomial, fe_oshape::OrientedShape, side_face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a finite element relative monomial vs. a side relative monomial.
integral_fe_rel_x_side_rel_on_side(mon1::Monomial, mon2::Monomial, fe_oshape::OrientedShape, side_face::FERelFace, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

#
####################################################################


# The following have default implementations or are defined in terms of the required
# generic functions and usually won't need to be implemented for specific mesh types.

function integral_face_rel_on_face(c::R, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh)
  if c == zeroR
    zeroR
  else
    c * integral_face_rel_on_face(one_mon(mesh), fe_oshape, face, mesh)
  end
end

function integral_face_rel_on_face(p::Polynomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_face_rel_on_face(p.mons[i], fe_oshape, face, mesh)
  end
  sum
end

integral_face_rel_on_face(pmc::NomialOrConst, fe::FENum, face::FERelFace, mesh::AbstractMesh) =
  integral_face_rel_on_face(pmc, Mesh.oriented_shape_for_fe(fe, mesh), face, mesh)

function integral_global_x_face_rel_on_fe_face(f::Function, p::Polynomial, fe::FENum, face::FERelFace, mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_global_x_face_rel_on_fe_face(f, p.mons[i], fe, face, mesh)
  end
  sum
end

# Implementations could specialize this to avoid creating the product monomial.
integral_face_rel_x_face_rel_on_face(mon1::Monomial, mon2::Monomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh) =
  integral_face_rel_on_face(mon1 * mon2, fe_oshape, face, mesh)

function integral_face_rel_x_face_rel_on_face(mon::Monomial, p::Polynomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_face_rel_x_face_rel_on_face(mon, p.mons[i], fe_oshape, face, mesh)
  end
  sum
end

integral_face_rel_x_face_rel_on_face(p::Polynomial, mon::Monomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh) =
  integral_face_rel_x_face_rel_on_face(mon, p, fe_oshape, face, mesh)

function integral_face_rel_x_face_rel_on_face(p1::Polynomial, p2::Polynomial, fe_oshape::OrientedShape, face::FERelFace, mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p1.coefs), j=1:length(p2.coefs)
    sum += p1.coefs[i] * p2.coefs[j] * integral_face_rel_x_face_rel_on_face(p1.mons[i], p2.mons[j], fe_oshape, face, mesh)
  end
  sum
end


# origin (see "Interpreting Monomials..." discussion at top).
function integral_side_rel_x_fe_rel_vs_outward_normal_on_side(p::Polynomial, q::VectorMonomial, fe_oshape::OrientedShape, side_face::FERelFace, mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.mons)
    sum += p.coefs[i] * integral_side_rel_x_fe_rel_vs_outward_normal_on_side(p.mons[i], q, fe_oshape, side_face, mesh)
  end
  sum
end

function fe_interior_origin(fe::FENum, mesh::AbstractMesh)
  const d = mesh.space_dim
  const coords = Array(R, d)
  fe_interior_origin!(fe, coords, mesh)
  coords
end

# TODO: unit tests
function num_non_boundary_sides_for_fe(fe::FENum, mesh::AbstractMesh)
  nb_sides = 0
  for sf=fe_face(1):Mesh.num_side_faces_for_fe(fe, mesh)
    if !Mesh.is_boundary_side(fe, sf, mesh) nb_sides += 1 end
  end
  nb_sides
end

function max_num_shape_sides(mesh::AbstractMesh)
  max_sides= 0
  for os=oshape(1):num_oriented_element_shapes(mesh)
    max_sides = max(max_sides, num_side_faces_for_shape(os, mesh))
  end
  max_sides
end

import Base.isequal
isequal(incls1::NBSideInclusions, incls2::NBSideInclusions) =
  incls1.fe1 == incls2.fe1 && incls1.face_in_fe1 == incls2.face_in_fe1 &&
  incls1.fe2 == incls2.fe2 && incls1.face_in_fe2 == incls2.face_in_fe2
import Base.hash
hash(incls::NBSideInclusions) =
  fe1 + 3*face_in_fe1 + 5 * fe2 + 7*face_in_fe2

end # end of module
