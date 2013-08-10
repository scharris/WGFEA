using Base.Test
using TMesh
using Common
import Mesh.NBSideInclusions, Mesh.FENum, Mesh.FEFaceNum, Mesh.oshapenum, Mesh.nbsidenum, Mesh.fenum, Mesh.fefacenum
import RMesh, RMesh.RectMesh
import Poly, Poly.Monomial, Poly.VectorMonomial


x = Monomial(1,0)
y = Monomial(0,1)
onemon = Monomial(0,0)

# Local monomials extended as global functions.

function local_mon_on_fe_int(exp1::Deg, exp2::Deg, fe::FENum, mesh::TriMesh)
  const o = Mesh.fe_interior_origin(fe, mesh)
  (x::Array{R,1}) -> (x[1]-o[1])^exp1 * (x[2]-o[2])^exp2
end

function local_mon_on_fe_side(exp1::Deg, exp2::Deg, fe::FENum, side_face::FEFaceNum, mesh::TriMesh)
  const ref_tri = TMesh.ref_tri_for_fe(fe, mesh)
  const fe_v1 = mesh.fes[fe].v1
  const a,b = TMesh.sideface_endpoint_pair(side_face,
                                           fe_v1, fe_v1 + ref_tri.v12, fe_v1 + ref_tri.v13,
                                           ref_tri.nums_faces_between_vertexes,
                                           false) # lesser endpoints first => false
  const o = 1/2*(a + b) # The local origin of each side face is the midpoint.
  (x::Array{R,1}) -> (x[1]-o._1)^exp1 * (x[2]-o._2)^exp2
end


# Test retrieval of side face endpoints.

# Return all side face endpoint pairs, with endpoints in traversal order.
@test TMesh.sideface_endpoint_pairs(Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (2,1,2), false) ==
  [(Point(0.,0.),Point(.5,.5)),
   (Point(.5,.5),Point(1.,1.)),
   (Point(1.,1.),Point(-1.,1.)),
   (Point(-1.,1.),Point(-.5,.5)),
   (Point(-.5,.5),Point(0.,0.))]

# Same as above but order endpoints within endpoint pairs in a standardized order with lesser point first.
@test TMesh.sideface_endpoint_pairs(Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (2,1,2), true) ==
  [(Point(0.,0.),Point(.5,.5)),
   (Point(.5,.5),Point(1.,1.)),
   (Point(-1.,1.), Point(1.,1.)),
   (Point(-1.,1.),Point(-.5,.5)),
   (Point(-.5,.5),Point(0.,0.))]

# Test individual side endpoint pair function, endpoints in traversal order.
@test [TMesh.sideface_endpoint_pair(fefacenum(sf), Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (2,1,2), false) for sf in 1:5] ==
  [(Point(0.,0.),Point(.5,.5)),
   (Point(.5,.5),Point(1.,1.)),
   (Point(1.,1.),Point(-1.,1.)),
   (Point(-1.,1.),Point(-.5,.5)),
   (Point(-.5,.5),Point(0.,0.))]

# Same as above but with endpoints in standardized order.
@test [TMesh.sideface_endpoint_pair(fefacenum(sf), Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (2,1,2), true) for sf in 1:5] ==
  [(Point(0.,0.),Point(.5,.5)),
   (Point(.5,.5),Point(1.,1.)),
   (Point(-1.,1.), Point(1.,1.)),
   (Point(-1.,1.),Point(-.5,.5)),
   (Point(-.5,.5),Point(0.,0.))]

# Check outward normals with number of faces between vertexes varying.
onormals = TMesh.outward_normals(Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (1,2,3))
@test onormals[1] == Vec(1/sqrt(2),-1/sqrt(2))
@test onormals[2] == Vec(0.,1.)
@test onormals[3] == Vec(0.,1.)
@test onormals[4] == Vec(-1/sqrt(2),-1/sqrt(2))
@test onormals[5] == Vec(-1/sqrt(2),-1/sqrt(2))
@test onormals[6] == Vec(-1/sqrt(2),-1/sqrt(2))


# Test integration between upper and lower lines converging to a point.

# Integrate the one monomial between lines of slope -3 and 2, diverging from point (-1,-2) to x = 3.
# horizontal span of generated triangle is 4
# according to the slopes, the right face of the resulting triangle has length 3*4 + 2*4 = 20
# area of triangle = 1/2 * 20 * 4 = 40
@test isapprox(TMesh.integral_mon_between_lines_meeting_at_point_and_vert_line(onemon, Point(-1.,-2.), -3.,2., 3.), 40., atol=1e-15)

# Monomial xy between lines diverging from point (2,1) of slope +/- 1, from x = 2 to x = 3.
# int_{x=2..3} int_{y=1+(x-2)..y=1-(x-2)} x y dy dx
#  = int_{x=2..3} [1/2 x y^2]_|{y=1+(x-2)..1-(x-2)}
#  = int_{x=2..3} x/2 ( (1+(x-2))^2 - (1-(x-2))^2 )  # (1+(x-2) + 1-(x-2)) (1+(x-2) - 1-(x-2)) = 2 (2(x-2)) = 4(x-2)
#  = int_{x=2..3} x/2 4 (x-2)
#  = int_{x=2..3} 2 x (x-2)
#  = int_{x=2..3} 2x^2 - 4x
#  = [2/3 x^3 - 2x^2]_|{x=2..3}
#  = 2 + 2/3
@test isapprox(TMesh.integral_mon_between_lines_meeting_at_point_and_vert_line(x*y, Point(2.,1.), -1.,1., 3.), 2 + 2/3, atol=1e-15)

# mirror image of the above about the y axis
@test isapprox(TMesh.integral_mon_between_lines_meeting_at_point_and_vert_line(x*y, Point(-2.,1.), -1.,1., -3.), -(2 + 2/3), atol=1e-15)

# input mesh, before subdivision
#     (0,3)
#       | .
#       |    .
#  (0,0)|_______.(4,0)
#
right_tri_mesh =
    """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    3
    1 0 0 0
    2 4 0 0
    3 0 3 0
    \$EndNodes
    \$Elements
    7
    1 15 2 0 1 1
    2 15 2 0 2 2
    3 15 2 0 3 3
    4 1 2 0 1 1 2
    5 1 2 0 2 2 3
    6 1 2 0 3 3 1
    7 2 2 0 5 1 2 3
    \$EndElements
    """

# one subdivison of input mesh
tmsh = TriMesh(IOString(right_tri_mesh), 1)

@test Mesh.space_dim(tmsh) === dim(2)
@test Mesh.one_mon(tmsh) == Poly.Monomial(0,0)
@test Mesh.num_fes(tmsh) === fenum(4)
@test Mesh.num_nb_sides(tmsh) === nbsidenum(3)
@test Mesh.num_oriented_element_shapes(tmsh) === oshapenum(2)
@test Mesh.oriented_shape_for_fe(fenum(1), tmsh) === oshapenum(1)
@test Mesh.num_side_faces_for_fe(fenum(1), tmsh) === fefacenum(3)
@test Mesh.num_side_faces_for_shape(oshapenum(1), tmsh) === fefacenum(3)
@test Mesh.dependent_dim_for_oshape_side(oshapenum(1), fefacenum(1), tmsh) === dim(2)

# 3 primary subtriangles
@test tmsh.fes[1] == ElTri(oshapenum(1),Vec(0.0,0.0))
@test let origin = [-1.,-1.]; Mesh.fe_interior_origin!(fenum(1), origin, tmsh); origin == [0.,0.] end
@test Mesh.fe_interior_origin(fenum(1), tmsh) == [0.,0.]
@test tmsh.fes[2] == ElTri(oshapenum(1),Vec(4/2,0.0))
@test Mesh.fe_interior_origin(fenum(2), tmsh) == [2.,0.]
@test tmsh.fes[3] == ElTri(oshapenum(1),Vec(0.0,3/2))
@test Mesh.fe_interior_origin(fenum(3), tmsh) == [0.,3/2]
# secondary subtriangle
@test tmsh.fes[4] == ElTri(oshapenum(2),Vec(4/2,0.0))
@test Mesh.fe_interior_origin(fenum(4), tmsh) == [2.,0.]
# primary reference triangle
@test tmsh.oshapes[1].v12 == Vec(2.,0.)
@test tmsh.oshapes[1].v13 == Vec(0.,1.5)
@test tmsh.oshapes[1].outward_normals_by_sideface == [Vec(0.,-1.), Vec(3/5,4/5), Vec(-1.,-0.)]
@test Mesh.shape_diameter_inv(oshapenum(1), tmsh) == 1/2.5
@test Mesh.dependent_dim_for_oshape_side(oshapenum(1), fefacenum(1), tmsh) === dim(2)
@test Mesh.dependent_dim_for_oshape_side(oshapenum(1), fefacenum(2), tmsh) === dim(2)
@test Mesh.dependent_dim_for_oshape_side(oshapenum(1), fefacenum(3), tmsh) === dim(1)
# secondary reference triangle
@test tmsh.oshapes[2].v12 == Vec(0.,1.5)
@test tmsh.oshapes[2].v13 == Vec(-2.,1.5)
@test tmsh.oshapes[2].outward_normals_by_sideface == [Vec(1.,-0.), Vec(0.,1.), Vec(-3/5,-4/5)]
@test Mesh.shape_diameter_inv(oshapenum(2), tmsh) == 1/2.5
@test Mesh.dependent_dim_for_oshape_side(oshapenum(2), fefacenum(1), tmsh) === dim(1)
@test Mesh.dependent_dim_for_oshape_side(oshapenum(2), fefacenum(2), tmsh) === dim(2)
@test Mesh.dependent_dim_for_oshape_side(oshapenum(2), fefacenum(3), tmsh) === dim(2)

@test Mesh.max_fe_diameter(tmsh) == 2.5

@test  Mesh.num_boundary_sides(tmsh) == 6
@test  Mesh.is_boundary_side(fenum(1), fefacenum(1), tmsh)
@test !Mesh.is_boundary_side(fenum(1), fefacenum(2), tmsh)
@test !Mesh.is_boundary_side(fenum(1), fefacenum(2), tmsh)
@test  Mesh.is_boundary_side(fenum(1), fefacenum(3), tmsh)
@test  Mesh.is_boundary_side(fenum(2), fefacenum(1), tmsh)
@test  Mesh.is_boundary_side(fenum(2), fefacenum(2), tmsh)
@test !Mesh.is_boundary_side(fenum(2), fefacenum(3), tmsh)
@test !Mesh.is_boundary_side(fenum(3), fefacenum(1), tmsh)
@test  Mesh.is_boundary_side(fenum(3), fefacenum(2), tmsh)
@test  Mesh.is_boundary_side(fenum(3), fefacenum(3), tmsh)


# Check non-boundary side inclusions in finite elements.
@test length(tmsh.nbsideincls_by_nbsidenum) == 3
fe1face2_nbsidenum = Mesh.nb_side_num_for_fe_side(fenum(1),fefacenum(2), tmsh)
@test tmsh.nbsidenums_by_feface[fenum(4),fefacenum(3)] == fe1face2_nbsidenum
@test Mesh.fe_inclusions_of_nb_side(fe1face2_nbsidenum, tmsh) ==
      NBSideInclusions(fe1face2_nbsidenum, fenum(1), fefacenum(2), fenum(4), fefacenum(3))
fe2face3_nbsidenum = Mesh.nb_side_num_for_fe_side(fenum(2),fefacenum(3), tmsh)
@test tmsh.nbsidenums_by_feface[fenum(4),fefacenum(1)] == fe2face3_nbsidenum
@test Mesh.fe_inclusions_of_nb_side(fe2face3_nbsidenum, tmsh) ==
      NBSideInclusions(fe2face3_nbsidenum, fenum(2), fefacenum(3), fenum(4), fefacenum(1))
fe3face1_nbsidenum = Mesh.nb_side_num_for_fe_side(fenum(3),fefacenum(1), tmsh)
@test tmsh.nbsidenums_by_feface[fenum(4),fefacenum(2)] == fe3face1_nbsidenum
@test Mesh.fe_inclusions_of_nb_side(fe3face1_nbsidenum, tmsh) ==
      NBSideInclusions(fe3face1_nbsidenum, fenum(3), fefacenum(1), fenum(4), fefacenum(2))


# Test integrals.

# primary reference element
#     (0,1.5)
#       | .
#       |    .
#  (0,0)|_______.(2,0)
#

# integrals of face-relative monomials on side faces

# integrals of constant 1 monomials should give side lengths
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(1), tmsh),
               2., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(2), tmsh),
               hypot(2.,1.5), rtol=1e-15, atol=1e-12)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(3), tmsh),
               1.5, rtol=1e-15, atol=1e-12)

@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)

# These compute side lengths like the above, but are expressed as a function times a local monomial.
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(1), tmsh),
               2., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(2), tmsh),
               hypot(2.,1.5), atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(3), tmsh),
               1.5, atol=1e-15, rtol=1e-15)

# integral of xy on side 1
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
xy_on_fe1_side1 = local_mon_on_fe_side(deg(1), deg(1), fenum(1), fefacenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side1, onemon, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)
x_on_fe1_side1 = local_mon_on_fe_side(deg(1), deg(0), fenum(1), fefacenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side1, y, fenum(1), fefacenum(1), tmsh),
               0., atol=1e-15, rtol=1e-15)

# Integrate polynomial xy on side face 2 of fe/oshape 1.
# We can traverse side 2 of fe 1's reference triangle in side relative coordinates with
# p(t) = (2,0)+t((0,1.5)-(2,0))-(1,0.75), 0<=t<=1. Then
# xy(p(t)) = xy((1 - 2t, 1.5t - 0.75)) = (1-2t)(1.5t-0.75)
# so the integral should be
# int_{0,..1} (1-2t)(1.5t-0.75) |(-2,1.5)| dt = 2.5[-.75t + 1.5t^2 - t^3]|{t=0,1} = -0.625
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
x_on_fe1_side2 = local_mon_on_fe_side(deg(1), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, onemon, fenum(1), fefacenum(2), tmsh),
               -0.625, atol=1e-15, rtol=1e-15)

# Integrate the same monomial over the same side but now also vs. outward normal.
one_comp1_vmon = VectorMonomial(onemon, dim(1))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y, one_comp1_vmon, oshapenum(1), fefacenum(2), tmsh),
               -0.625 * 3/5, atol=1e-15, rtol=1e-15) # 3/5 = component 1 of outward normal
one_comp2_vmon = VectorMonomial(onemon, dim(2))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x*y, one_comp2_vmon, oshapenum(1), fefacenum(2), tmsh),
               -0.625 * 4/5, atol=1e-15, rtol=1e-15) # 4/5 = component 2 of outward normalj
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, one_comp1_vmon, oshapenum(1), fefacenum(1), tmsh),
               0, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, one_comp2_vmon, oshapenum(1), fefacenum(3), tmsh),
               0, atol=1e-15, rtol=1e-15)

# Integrate (0,x) vector monomial vs outward normal along first side, which should equal int_0^2 -x dx = -2.
x_comp2_vmon = VectorMonomial(x, dim(2))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(onemon, x_comp2_vmon, oshapenum(1), fefacenum(1), tmsh),
               -2, atol=1e-15, rtol=1e-15)

# Like the above but with a scalar monomial instead of vector dotted against the normal.
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x, onemon, oshapenum(1), fefacenum(1), tmsh),
               2, atol=1e-15, rtol=1e-15)

# Integrate x side-local monomial vs (0,x) vector monomial (fe-relative) vs outward normal along the first side.
# As an fe-local monomial, the side local monomial x is x-1. Thus the integral should be
# int_0^2 -(x-1)x dx = -1/3 2^3 + 2^2/2 = -2/3
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(x, x_comp2_vmon, oshapenum(1), fefacenum(1), tmsh),
               -2/3, atol=1e-15, rtol=1e-15)
# Like the above but with a scalar monomial instead of vector dotted against the normal.
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(x, x, oshapenum(1), fefacenum(1), tmsh),
               2/3, atol=1e-15, rtol=1e-15)

# Integrate y side-local monomial vs (y,0) fe-relative vector monomial vs outward normal along the third side.
# As an fe-local monomial, the side local monomial y is y-3/4. Thus the integral should be
# int_0^{3/2} -(y-3/4)y dy = -1/2 (3/2)^2 + 3/8 (3/2)^2 = -9/32
y_comp1_vmon = VectorMonomial(y, dim(1))
@test isapprox(Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(y, y_comp1_vmon, oshapenum(1), fefacenum(3), tmsh),
               -9/32, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_fe_rel_x_side_rel_on_oshape_side(y, y, oshapenum(1), fefacenum(3), tmsh),
               9/32, atol=1e-15, rtol=1e-15)


# interior integrals

# Integral of constant one monomial on a triangle interior should give its area.
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*2*1.5, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*2*1.5, atol=1e-15, rtol=1e-15)

# integral xy^2 over the interior of reference triangle 1
# = int_{x=0..2} int_{y=0..1.5-0.75x} xy^2 = int_{x=0..2} 1/3 x(1.5-0.75x)^3
# = 1/3 (1.6875 x^2 + -1.6875 x^3 + 0.6328125 x^4 + -0.084375 x^5)|_{x=0,2}
# = 0.225
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y^2, oshapenum(1), fefacenum(0), tmsh),
               0.225, atol=1e-15, rtol=1e-15)


# Integrate function (x,y)->(x,y)-o(fe) against local monomial y on various finite element interiors, which
# should equal the integral of xy^2 local monomial over the oriented shape's interior.
xy_on_fe1 = local_mon_on_fe_int(deg(1),deg(1),fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)
xy_on_fe2 = local_mon_on_fe_int(deg(1),deg(1),fenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe2, y, fenum(2), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)
xy_on_fe3 = local_mon_on_fe_int(deg(1),deg(1),fenum(3), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe3, y, fenum(3), fefacenum(0), tmsh),
               0.225, rtol=1e-13, atol=1e-13)

# input mesh, before subdivision
#     (2,8).
#         . .
#        .   .
#  (0,0). . . .(4,0)
#
isosc_tri_mesh =
    """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    3
    1 0 0 0
    2 4 0 0
    3 2 8 0
    \$EndNodes
    \$Elements
    7
    1 15 2 0 1 1
    2 15 2 0 2 2
    3 15 2 0 3 3
    4 1 2 0 1 1 2
    5 1 2 0 2 2 3
    6 1 2 0 3 3 1
    7 2 2 0 5 1 2 3
    \$EndElements
    """

# Check that integrals that must be done in two pieces are done properly.
tmsh = TriMesh(IOString(isosc_tri_mesh), 2)
@test tmsh.fes[1] == ElTri(oshapenum(1),Vec(0.,0.))

# primary reference element
#   (1/2,2).
#         . .
#        .   .
#  (0,0). . . .(1,0)
#

@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)

# Integrate the local monomial xy^2 over the interior of a reference triangle.
# int_0^2 int_{y/4}^{-y/4 + 1}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((-y/4 + 1)^2 - (y/4)^2) dy
#  = int_0^2 1/2 y^2 - 1/4 y^3 dy
#  = (1/6 y^3 - 1/16 y^4)|_0^2
#  = 8/6 - 1 = 1/3
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(1), fefacenum(0), tmsh),
               1/3, rtol=1e-11, atol=1e-11)
# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe1 = local_mon_on_fe_int(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               1/3, rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate xy^2 monomial along side face 1, on which y = 0.
xy_on_fe1 = local_mon_on_fe_int(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(1), tmsh),
               0., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y^2, fenum(1), fefacenum(1), tmsh),
               0., rtol=1e-11, atol=1e-11)

# Integrate polynomial xy on side face 2 of fe/oshape 1.
# We can traverse side 2 of fe 1's reference triangle in side relative coordinates with
# p(t) = (1,0) + t((1/2,2)-(1,0)) - (0.75,1)
#      = (1/4,-1) + t(-1/2,2),   (0<=t<=1). Then
# xy(p(t)) = xy((1/4 - 1/2 t, -1 + 2t)) = -1/4 + t - t^2
# so the integral should be
# int_0^1 (-1/4 + t - t^2) |(-1/2,2)| dt = sqrt(17)/2 [-t/4 + 1/2 t^2 -1/3 t^3]|_{t=0,1} = -sqrt(17)/24
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x*y, oshapenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
x_on_fe1_side2 = local_mon_on_fe_side(deg(1), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(x_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, onemon, fenum(1), fefacenum(2), tmsh),
               -sqrt(17)/24, atol=1e-15, rtol=1e-15)


# inversion of the above, original undivided triangle
#  (-2,8). . . .(2,8)
#         .   .
#          . .
#           .(0,0)
inv_isosc_tri_mesh =
    """
    \$MeshFormat
    2.2 0 8
    \$EndMeshFormat
    \$Nodes
    3
    1 0 0 0
    2 2 8 0
    3 -2 8 0
    \$EndNodes
    \$Elements
    7
    1 15 2 0 1 1
    2 15 2 0 2 2
    3 15 2 0 3 3
    4 1 2 0 1 1 2
    5 1 2 0 2 2 3
    6 1 2 0 3 3 1
    7 2 2 0 5 1 2 3
    \$EndElements
    """

# primary reference element
# (-1/2,2). . . .(1/2,2)
#          .   .
#           . .
#            .(0,0)

# Check that integrals that must be done in two pieces are done properly.
tmsh = TriMesh(IOString(inv_isosc_tri_mesh), 2)
@test tmsh.fes[1] == ElTri(oshapenum(1),Vec(0.,0.))

@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(1), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_face_rel_on_oshape_face(onemon, oshapenum(2), fefacenum(0), tmsh),
               0.5*1*2, atol=1e-15, rtol=1e-15)

# Integrate the local monomial xy^2 over the interior of a reference triangle.
# int_0^2 int_{-y/4}^{y/4}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((y/4)^2 - (-y/4)^2) dy
#  = 0
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(1), fefacenum(0), tmsh),
               0., rtol=1e-11, atol=1e-11)

# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe1 = local_mon_on_fe_int(deg(1),deg(1), fenum(1), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1, y, fenum(1), fefacenum(0), tmsh),
               0., rtol=1e-11, rtol=1e-11)

@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate the local monomial x^2y^2 over the interior of the primary reference triangle.
# int_0^2 int_{-y/4}^{y/4}  x^2y^2 dx dy
#  = int_0^2 1/3 y^2 ((y/4)^3 - (-y/4)^3) dy
#  = int_0^2 1/3 y^2 2/64 y^3 dy
#  = int_0^2 1/96 y^5 dy
#  = 1/96 1/6 y^6|_0^2
#  = 1/9
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x^2 * y^2, oshapenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)

# Same as above, but expressed as a product of a global function and monomial.

@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x^2 * y^2, fenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(local_mon_on_fe_int(deg(2),deg(2), fenum(1), tmsh), onemon, fenum(1), fefacenum(0), tmsh),
               1/9, rtol=1e-11, atol=1e-11)

# areas
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(1), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(2), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(3), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., onemon, fenum(4), fefacenum(0), tmsh),
               1., rtol=1e-11, atol=1e-11)

# Integrate polynomial x^3y on side face 2 of fe/oshape 1.
# The integral is int_{-1/2}^{1/2} 2 x^3 dx = 0.
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x^3*y, oshapenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x^3*y, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
xcubed_on_fe1_side2 = local_mon_on_fe_side(deg(3), deg(0), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xcubed_on_fe1_side2, y, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)
xy_on_fe1_side2 = local_mon_on_fe_side(deg(1), deg(1), fenum(1), fefacenum(2), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe1_side2, x^2, fenum(1), fefacenum(2), tmsh),
               0, atol=1e-15, rtol=1e-15)


# secondary reference element
#  (-1/2,2).
#         . .
#        .   .
# (-1,0). . . .(0,0)

# Integrate the local monomial xy^2 over the interior of the secondary reference triangle.
# int_0^2 int_{y/4-1}^{-y/4}  xy^2 dx dy
#  = int_0^2 1/2 y^2 ((-y/4)^2 - (y/4-1)^2) dy
#  = int_0^2 1/4 y^3 - 1/2 y^2 dy
#  = (1/16 y^4 - 1/6 y^3)|_0^2
#  = 1 - 8/6 = -1/3
@test isapprox(Mesh.integral_face_rel_on_oshape_face(x * y^2, oshapenum(2), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)

# Like the above, but expressed as a product of a global function, chosen to match local monomial xy, and monomial y.
xy_on_fe4 = local_mon_on_fe_int(deg(1),deg(1), fenum(4), tmsh)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(xy_on_fe4, y, fenum(4), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)
@test isapprox(Mesh.integral_global_x_face_rel_on_fe_face(_->1., x*y^2, fenum(4), fefacenum(0), tmsh),
               -1/3, atol=1e-11, rtol=1e-11)


# Test construction from rectangle mesh.
rmesh2x2 = RectMesh([1.,1.], [3.,3.], [RMesh.mesh_coord(2), RMesh.mesh_coord(2)])
tmesh = from_rect_mesh(rmesh2x2, 1)
@test Mesh.num_fes(tmesh) == 32
@test Mesh.num_boundary_sides(tmesh) == 16
