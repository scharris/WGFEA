using Base.Test
using TMesh
using Common
import Mesh.NBSideInclusions, Mesh.oshapenum, Mesh.nbsidenum, Mesh.fenum, Mesh.fefacenum

one_tri_mesh =
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


tmsh0 = TriMesh(IOString(one_tri_mesh), 0)
@test tmsh0.fes == [ElTri(oshapenum(1), Vec(0.,0.))]
@test tmsh0.oshapes[1].v12 == Vec(4.,0.)
@test tmsh0.oshapes[1].v13 == Vec(0.,3.)
@test tmsh0.oshapes[1].nums_faces_between_vertexes == (1,1,1)
@test tmsh0.oshapes[1].outward_normals_by_sideface == [Vec(0.,-1.), Vec(3/5,4/5), Vec(-1.,-0.)]
@test tmsh0.oshapes[1].diameter_inv == 1/5

@test length(tmsh0.nbsidenums_by_feface) == 0
@test tmsh0.dep_dims_by_nbsidenum == []
@test tmsh0.num_b_sides == 3

tmsh1 = TriMesh(IOString(one_tri_mesh), 1)
# 3 primary subtriangles
@test tmsh1.fes[1] == ElTri(oshapenum(1),Vec(0.0,0.0))
@test tmsh1.fes[2] == ElTri(oshapenum(1),Vec(4/2,0.0))
@test tmsh1.fes[3] == ElTri(oshapenum(1),Vec(0.0,3/2))
# secondary subtriangle
@test tmsh1.fes[4] == ElTri(oshapenum(2),Vec(4/2,0.0))
# primary reference triangle
@test tmsh1.oshapes[1].v12 == Vec(2.,0.)
@test tmsh1.oshapes[1].v13 == Vec(0.,1.5)
@test tmsh1.oshapes[1].outward_normals_by_sideface == [Vec(0.,-1.), Vec(3/5,4/5), Vec(-1.,-0.)]
@test tmsh1.oshapes[1].diameter_inv == 1/2.5
# secondary reference triangle
@test tmsh1.oshapes[2].v12 == Vec(0.,1.5)
@test tmsh1.oshapes[2].v13 == Vec(-2.,1.5)
@test tmsh1.oshapes[2].outward_normals_by_sideface == [Vec(1.,-0.), Vec(0.,1.), Vec(-3/5,-4/5)]
@test tmsh1.oshapes[2].diameter_inv == 1/2.5

@test tmsh1.num_b_sides == 6

# Check non-boundary side inclusions in finite elements.
@test length(tmsh1.nbsideincls_by_nbsidenum) == 3
fe1face2_nbsidenum = tmsh1.nbsidenums_by_feface[(fenum(1),fefacenum(2))]
@test tmsh1.nbsidenums_by_feface[(fenum(4),fefacenum(3))] == fe1face2_nbsidenum
@test tmsh1.nbsideincls_by_nbsidenum[fe1face2_nbsidenum] ==
      NBSideInclusions(fe1face2_nbsidenum, fenum(1), fefacenum(2), fenum(4), fefacenum(3))
fe2face3_nbsidenum = tmsh1.nbsidenums_by_feface[(fenum(2),fefacenum(3))]
@test tmsh1.nbsidenums_by_feface[(fenum(4),fefacenum(1))] == fe2face3_nbsidenum
@test tmsh1.nbsideincls_by_nbsidenum[fe2face3_nbsidenum] ==
      NBSideInclusions(fe2face3_nbsidenum, fenum(2), fefacenum(3), fenum(4), fefacenum(1))
fe3face1_nbsidenum = tmsh1.nbsidenums_by_feface[(fenum(3),fefacenum(1))]
@test tmsh1.nbsidenums_by_feface[(fenum(4),fefacenum(2))] == fe3face1_nbsidenum
@test tmsh1.nbsideincls_by_nbsidenum[fe3face1_nbsidenum] ==
      NBSideInclusions(fe3face1_nbsidenum, fenum(3), fefacenum(1), fenum(4), fefacenum(2))

# Check dependent dimensions for non-boundary sides.
@test tmsh1.dep_dims_by_nbsidenum[fe1face2_nbsidenum] == dim(2)
@test tmsh1.dep_dims_by_nbsidenum[fe2face3_nbsidenum] == dim(1)
@test tmsh1.dep_dims_by_nbsidenum[fe3face1_nbsidenum] == dim(2)


# Check outward normals with number of faces between vertexes varying.
onormals = TMesh.outward_normals(Vec(0.,0.), Vec(1.,1.), Vec(-1.,1.), (1,2,3))
@test onormals[1] == Vec(1/sqrt(2),-1/sqrt(2))
@test onormals[2] == Vec(0.,1.)
@test onormals[3] == Vec(0.,1.)
@test onormals[4] == Vec(-1/sqrt(2),-1/sqrt(2))
@test onormals[5] == Vec(-1/sqrt(2),-1/sqrt(2))
@test onormals[6] == Vec(-1/sqrt(2),-1/sqrt(2))

