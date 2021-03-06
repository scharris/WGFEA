module Mesh

export FENum, fenum,
       NBSideNum, nbsidenum,
       FEFaceNum, fefacenum, interior_face, feface_one,
       OShapeNum, oshapenum,
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
       dependent_dim_for_oshape_side,
       fe_inclusions_of_nb_side,
       nb_side_num_for_fe_side,
       is_boundary_side,
       num_boundary_sides,
       num_non_boundary_sides_for_fe,
       integral_face_rel_on_oshape_face,
       integral_global_x_face_rel_on_fe_face,
       integral_face_rel_x_face_rel_on_oshape_face,
       integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side,
       integral_fe_rel_x_side_rel_on_oshape_side

require("Common.jl")
require("Poly.jl")

using Common
import Poly.Monomial, Poly.Polynomial, Poly.VectorMonomial

# Number type for enumeration of all finite elements in a mesh.
typealias FENum Uint64
fenum(i::Integer) = if i > typemax(FENum) || i <= 0 error("fe number out of range") else convert(FENum, i) end

# Number type for enumeration of all non-boundary sides in a mesh.
typealias NBSideNum Uint64
nbsidenum(i::Integer) = if i > typemax(NBSideNum) || i <= 0 error("nb side number out of range") else convert(NBSideNum, i) end

# Finite Element Relative Face - Number type for enumerating the interior and
# side parts of a finite element relative to the element itself. For all
# meshes, the interior face is always numbered 0, while the sides are numbered
# 1 ... num_side_faces_for_fe(fe, mesh) for a given fe.
typealias FEFaceNum Uint8
fefacenum(i::Integer) = if i >= typemax(FEFaceNum) || i < 0 error("fe face number out of range") else convert(FEFaceNum, i) end
const interior_face = fefacenum(0)
const feface_one = fefacenum(1) # first side face

typealias OShapeNum Uint16
oshapenum(i::Integer) = if i > typemax(OShapeNum) || i <= 0 error("oriented shape number out of range") else convert(OShapeNum, i) end
const oshape_one = oshapenum(1)

# Given a side not in the outside boundary, this structure represents the two
# finite elements which include the side, together with the side's face number
# in each of the finite elements.
immutable NBSideInclusions
  nb_side_num::NBSideNum
  fe1::FENum
  face_in_fe1::FEFaceNum
  fe2::FENum
  face_in_fe2::FEFaceNum
end


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

num_side_faces_for_fe(fe::FENum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

num_side_faces_for_shape(fe_oshape::OShapeNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")


# The dependent dim is, for a given side, the dimension j which is affine-
# dependent on the other dimensions on the side. That is, the function
# returns a j for which c_0,...,c_d exist such that
#             x_j = c_0 + sum_{i=1..d} c_i x_i for all (x_1,...,x_d) in the side.
# There may be more than one such coordinate number, in which case any one of
# these is returned.
dependent_dim_for_oshape_side(fe_oshape::OShapeNum, side_face::FEFaceNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

fe_inclusions_of_nb_side(side_num::NBSideNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Return non-boundary side number of the indicated fe relative side.
nb_side_num_for_fe_side(fe::FENum, side_face::FEFaceNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")


is_boundary_side(fe::FENum, face::FEFaceNum, mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

num_boundary_sides(mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

shape_diameter_inv(fe_oshape::OShapeNum, mesh::AbstractMesh) =
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
integral_face_rel_on_oshape_face(m::Monomial,
                                 fe_oshape::OShapeNum, face::FEFaceNum,
                                 mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a global function f on the indicated finite element face, multiplied against
# a monomial interpreted locally on the face.
integral_global_x_face_rel_on_fe_face(f::Function,
                                      mon::Monomial,
                                      fe::FENum, face::FEFaceNum,
                                      mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a side-local monomial vs. a vector monomial interpreted relative to the
# entire finite element, dot multiplied with the outward normal for the side.
integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(v::Monomial,
                                                            q::VectorMonomial,
                                                            fe_oshape::OShapeNum, side_face::FEFaceNum,
                                                            mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

# Integrate a finite element relative monomial vs. a side relative monomial.
integral_fe_rel_x_side_rel_on_oshape_side(mon1::Monomial,
                                          mon2::Monomial,
                                          fe_oshape::OShapeNum, side_face::FEFaceNum,
                                          mesh::AbstractMesh) =
  error("not implemented, mesh implementation is incomplete")

#
####################################################################


# The following have default implementations or are defined in terms of the required
# generic functions and usually won't need to be implemented for specific mesh types.

function integral_face_rel_on_oshape_face(c::R,
                                          fe_oshape::OShapeNum, face::FEFaceNum,
                                          mesh::AbstractMesh)
  if c == zeroR
    zeroR
  else
    c * integral_face_rel_on_oshape_face(one_mon(mesh), fe_oshape, face, mesh)
  end
end

function integral_face_rel_on_oshape_face(p::Polynomial,
                                          fe_oshape::OShapeNum, face::FEFaceNum,
                                          mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_face_rel_on_oshape_face(p.mons[i], fe_oshape, face, mesh)
  end
  sum
end


function integral_global_x_face_rel_on_fe_face(f::Function,
                                               p::Polynomial,
                                               fe::FENum, face::FEFaceNum,
                                               mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_global_x_face_rel_on_fe_face(f, p.mons[i], fe, face, mesh)
  end
  sum
end

# Implementations could specialize this to avoid creating the product monomial.
integral_face_rel_x_face_rel_on_oshape_face(mon1::Monomial,
                                            mon2::Monomial,
                                            fe_oshape::OShapeNum, face::FEFaceNum,
                                            mesh::AbstractMesh) =
  integral_face_rel_on_oshape_face(mon1 * mon2, fe_oshape, face, mesh)

function integral_face_rel_x_face_rel_on_oshape_face(mon::Monomial,
                                                     p::Polynomial,
                                                     fe_oshape::OShapeNum, face::FEFaceNum,
                                                     mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.coefs)
    sum += p.coefs[i] * integral_face_rel_x_face_rel_on_oshape_face(mon, p.mons[i], fe_oshape, face, mesh)
  end
  sum
end

integral_face_rel_x_face_rel_on_oshape_face(p::Polynomial,
                                            mon::Monomial,
                                            fe_oshape::OShapeNum, face::FEFaceNum,
                                            mesh::AbstractMesh) =
  integral_face_rel_x_face_rel_on_oshape_face(mon, p, fe_oshape, face, mesh)

function integral_face_rel_x_face_rel_on_oshape_face(p1::Polynomial,
                                                     p2::Polynomial,
                                                     fe_oshape::OShapeNum, face::FEFaceNum,
                                                     mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p1.coefs), j=1:length(p2.coefs)
    sum += p1.coefs[i] * p2.coefs[j] *
           integral_face_rel_x_face_rel_on_oshape_face(p1.mons[i], p2.mons[j], fe_oshape, face, mesh)
  end
  sum
end


# origin (see "Interpreting Monomials..." discussion at top).
function integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(p::Polynomial,
                                                                     q::VectorMonomial,
                                                                     fe_oshape::OShapeNum, side_face::FEFaceNum,
                                                                     mesh::AbstractMesh)
  sum = zeroR
  for i=1:length(p.mons)
    sum += p.coefs[i] *
           integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(p.mons[i], q, fe_oshape, side_face, mesh)
  end
  sum
end

function fe_interior_origin(fe::FENum, mesh::AbstractMesh)
  const d = space_dim(mesh)
  const coords = Array(R, d)
  fe_interior_origin!(fe, coords, mesh)
  coords
end

# TODO: unit tests
function num_non_boundary_sides_for_fe(fe::FENum, mesh::AbstractMesh)
  nb_sides = 0
  for sf=feface_one:num_side_faces_for_fe(fe, mesh)
    if !is_boundary_side(fe, sf, mesh) nb_sides += 1 end
  end
  nb_sides
end

function max_num_shape_sides(mesh::AbstractMesh)
  max_sides= 0
  for os=oshape_one:num_oriented_element_shapes(mesh)
    max_sides = max(max_sides, num_side_faces_for_shape(os, mesh))
  end
  max_sides
end

end # end of module
