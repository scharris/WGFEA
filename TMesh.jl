module TMesh
export TMesh

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature

typealias Point (R,R)
typealias PointNum Uint32

immutable RefTri
  # implicit first vertex at (0,0)
  vert_2::Point
  vert_3::Point
  side_dep_dims::Array{Dim, 1}
end

immutable ElTri
  oshape::OShapeNum
  vert_1::Point
end


end # end of module