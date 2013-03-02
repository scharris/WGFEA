module RMesh

export FENum, fe_num, no_fe,
       FECoord, fe_coord,
       FaceNum, face_num, interior_face, top_face, right_face, bottom_face, left_face,
       RectMesh,
       fe_row, fe_col,
       fe_left_of, fe_right_of, fe_above, fe_below

import Poly.R

# finite element number type
typealias FENum Uint64
fe_num(i::Integer) = uint64(i)
const no_fe = fe_num(0)

# finite element number type
typealias FECoord Uint32
fe_coord(i::Integer) = uint32(i)

typealias FaceNum Uint8
face_num(i::Integer) = if 0 <= i <= 4 uint8(i) else "face number out of range" end
const interior_face = face_num(0)
const top_face      = face_num(1)
const right_face    = face_num(2)
const bottom_face   = face_num(3)
const left_face     = face_num(4)
const no_face       = 0xff

type RectMesh
  bottom_left::(R,R)
  top_right::(R,R)
  rows::FECoord
  cols::FECoord
  # Computed values
  fe_width::R
  fe_height::R

  function RectMesh(bottom_left::(R,R), top_right::(R,R), nrows::FECoord, ncols::FECoord)
    fe_w = (top_right[1] - bottom_left[1]) / ncols
    fe_h = (top_right[2] - bottom_left[2]) / nrows

    new(bottom_left, top_right, nrows, ncols,
        fe_w, fe_h)
  end # RectMesh constructor

end # type RectMesh


function ==(rmesh1::RectMesh, rmesh2::RectMesh)
  rmesh1.bottom_left == rmesh2.bottom_left &&
  rmesh1.top_right   == rmesh2.top_right &&
  rmesh1.rows        == rmesh2.rows &&
  rmesh1.cols        == rmesh2.cols
end


# Returns the number of the finite element to the left, or 0 if there is none.
fe_left_of(fe::FENum) =
  let col = fe_col(fe)
    if col == 1 0 else fe - 1 end
  end

fe_row(fe::FENum) = div(fe-1, ncols) + 1
fe_col(fe::FENum) = mod(fe-1, ncols) + 1

# Returns the number of the finite element to the right, or 0 if there is none.
fe_right_of(fe::FENum) =
  let col = fe_col(fe)
    if col == ncols 0 else fe + 1 end
  end

# Returns the number of the finite element above the passed finite element, or 0 if there is none.
fe_above(fe::FENum) =
  let row = fe_row(fe)
    if row == nrows 0 else fe + ncols end
  end

# Returns the number of the finite element below the passed finite element, or 0 if there is none.
fe_below(fe::FENum) =
  let row = fe_row(fe)
    if row == 1 0 else fe - ncols end
  end

end # end of module
