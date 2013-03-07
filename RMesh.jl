module RMesh

export FENum, fe_num, no_fe,
       MeshRC, mesh_rc,
       FaceNum, face_num, interior_face, top_face, right_face, bottom_face, left_face,
       RectMesh,
       fe_row, fe_col,
       fe_rix, fe_cix,
       fe_coords!, fe_coords,
       fe_left_of, fe_right_of, fe_above, fe_below

import Poly.R

# finite element number type
typealias FENum Uint64
fe_num(i::Integer) = uint64(i)
const no_fe = fe_num(0)

# type for mesh row and column numbers
typealias MeshRC Uint32
mesh_rc(i::Integer) = uint32(i)

typealias FaceNum Uint8
face_num(i::Integer) = if 0 <= i <= 4 uint8(i) else "face number out of range" end
const interior_face = face_num(0)
const top_face      = face_num(1)
const right_face    = face_num(2)
const bottom_face   = face_num(3)
const left_face     = face_num(4)
const no_face       = 0xff

type RectMesh
  bottom_left_x::R
  bottom_left_y::R
  top_right_x::R
  top_right_y::R
  rows::MeshRC
  cols::MeshRC
  # Computed values
  fe_width::R
  fe_height::R

  function RectMesh(bottom_left::(R,R), top_right::(R,R), nrows::MeshRC, ncols::MeshRC)
    fe_w = (top_right[1] - bottom_left[1]) / ncols
    fe_h = (top_right[2] - bottom_left[2]) / nrows

    new(bottom_left[1], bottom_left[2], top_right[1], top_right[2], nrows, ncols,
        fe_w, fe_h)
  end # RectMesh constructor

end # type RectMesh


function ==(mesh1::RectMesh, mesh2::RectMesh)
  mesh1.bottom_left_x == mesh2.bottom_left_x &&
  mesh1.bottom_left_y == mesh2.bottom_left_y &&
  mesh1.top_right_x   == mesh2.top_right_x &&
  mesh1.top_right_y   == mesh2.top_right_y &&
  mesh1.rows          == mesh2.rows &&
  mesh1.cols          == mesh2.cols
end


fe_row(fe::FENum, mesh::RectMesh) = fe_rix(fe-1, mesh) + 1
fe_col(fe::FENum, mesh::RectMesh) = fe_cix(fe-1, mesh) + 1

fe_rix(fe_ix::FENum, mesh::RectMesh) = div(fe_ix, mesh.cols)
fe_cix(fe_ix::FENum, mesh::RectMesh) = mod(fe_ix, mesh.cols)

# Fill the passed 4 element array with the finite element rectangle coordingates,
# in the order: [bottom_left_x, bottom_left_y, top_right_x, top_right_y]
function fe_coords!(fe::FENum, mesh::RectMesh, bl_tr_coords::Vector{R})
  fe_ix = fe - 1
  fe_w = mesh.fe_width
  fe_h = mesh.fe_height
  rix = fe_rix(fe_ix, mesh)
  cix = fe_cix(fe_ix, mesh)
  bl_x = mesh.bottom_left_x + cix * fe_w
  bl_y = mesh.bottom_left_y + rix * fe_h
  tr_x = bl_x + fe_w
  tr_y = bl_y + fe_h
  bl_tr_coords[1] = bl_x
  bl_tr_coords[2] = bl_y
  bl_tr_coords[3] = tr_x
  bl_tr_coords[4] = tr_y
end

# Functional variant of the above.
fe_coords(fe::FENum, mesh::RectMesh) =
  let a = zeros(R,4)
    fe_coords!(fe, mesh, a)
    a
  end


# Returns the number of the finite element to the left, or 0 if there is none.
fe_left_of(fe::FENum, mesh::RectMesh) =
  let col = fe_col(fe, mesh)
    if col == 1 0 else fe - 1 end
  end

# Returns the number of the finite element to the right, or 0 if there is none.
fe_right_of(fe::FENum, mesh::RectMesh) =
  let col = fe_col(fe, mesh)
    if col == mesh.cols 0 else fe + 1 end
  end

# Returns the number of the finite element above the passed finite element, or 0 if there is none.
fe_above(fe::FENum, mesh::RectMesh) =
  let row = fe_row(fe, mesh)
    if row == mesh.rows 0 else fe + mesh.cols end
  end

# Returns the number of the finite element below the passed finite element, or 0 if there is none.
fe_below(fe::FENum, mesh::RectMesh) =
  let row = fe_row(fe, mesh)
    if row == 1 0 else fe - mesh.cols end
  end

end # end of module
