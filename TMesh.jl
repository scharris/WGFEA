module TMesh
export TMesh

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.RelFace, Mesh.OrientedShape, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.rface
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

# TODO
typealias PointNum Uint32

immutable Element
  points::(PointNum,PointNum,Pointnum)
  oshape::OrientedShape
end

immutable NBSide
  lesserEl::Element
  relFaceInLesserEl::RelFace
  greaterEl::Element
  relFaceInGreaterEl::RelFace
end

immutable TriangleMesh
  points::Array{(R,R),1}
  elements::Array{Element,1}
  nb_sides::Array{NBSide,1}
  b_sides::Array{BSide,1}
  # (fe,rface) -> 
end