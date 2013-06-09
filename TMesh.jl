module TMesh
export TriMesh,
       Vec,
       Point,
       RefTri,
       ElTri

using Common
import Mesh, Mesh.FENum, Mesh.NBSideNum, Mesh.FEFaceNum, Mesh.OShapeNum, Mesh.AbstractMesh,
       Mesh.NBSideInclusions, Mesh.fefacenum, Mesh.fenum, Mesh.nbsidenum
import Poly, Poly.Monomial, Poly.VectorMonomial
import Cubature.hcubature


# vector type for representing relative node offsets and points in the mesh
immutable Vec
  _1::R
  _2::R
end
typealias Point Vec

# reference triangle
immutable RefTri
  v12::Vec
  v13::Vec
  nums_faces_between_vertexes::(Int,Int,Int) # v1 v2, v2 v3, v3 v1
  # derived data
  outward_normals_by_sideface::Array{Vec,1}
  dep_dims_by_sideface::Array{Dim,1}
  diameter_inv::R
  num_side_faces::FEFaceNum

  function RefTri(v1::Vec, v2::Vec, v3::Vec, scale::R, nums_faces_btw_verts::(Int,Int,Int))
    const scaled_v12 = scale*(v2-v1)
    const scaled_v13 = scale*(v3-v1)
    const norm_scaled_v23 = scale * norm(v3-v2)
    const normals = outward_normals(v1,v2,v3, nums_faces_btw_verts)
    const dep_dims_by_sideface = sideface_dep_dims(v1,v2,v3, nums_faces_btw_verts)
    const diameter_inv = 1./max(norm(scaled_v12), norm(scaled_v13), norm_scaled_v23)
    const num_side_faces = sum(nums_faces_btw_verts)
    new(scaled_v12,
        scaled_v13,
        nums_faces_btw_verts,
        normals,
        dep_dims_by_sideface,
        diameter_inv,
        num_side_faces)
  end
end

# finite element triangle
immutable ElTri
  oshapenum::OShapeNum
  v1::Point
end

immutable TriMesh <: AbstractMesh
  # finite elements, indexed by finite element number
  fes::Array{ElTri,1}

  # reference triangles, indexed by oriented shape number
  oshapes::Array{RefTri,1}

  # non-boundary side numbers by fe face.
  nbsidenums_by_feface::Dict{(FENum,FEFaceNum), NBSideNum}

  # inclusions of non-boundary sides in fe's, indexed by non-boundary side number
  nbsideincls_by_nbsidenum::Array{NBSideInclusions,1}

  # number of boundary sides
  num_b_sides::Int64

  # cached items / computations / conversions
  space_dim::Dim
  one_mon::Monomial
  num_fes::FENum
  num_nb_sides::NBSideNum
  num_oshapes::OShapeNum

  # integration support members
  integrand_work_array::Vector{R}
  space_dim_zeros::Vector{R}
  space_dim_ones::Vector{R}
  integration_rel_err::R
  integration_abs_err::R
end

const zeroPt = Point(0.,0.)
const singleZero = [zeroR]
const singleOne = [oneR]

####################################################################
# implementations of required AbstractMesh functions
import Mesh.space_dim
space_dim(mesh::TriMesh) = mesh.space_dim

import Mesh.one_mon
one_mon(mesh::TriMesh) = mesh.one_mon

import Mesh.num_fes
num_fes(mesh::TriMesh) = mesh.num_fes

import Mesh.num_nb_sides
num_nb_sides(mesh::TriMesh) = mesh.num_nb_sides

import Mesh.num_oriented_element_shapes
num_oriented_element_shapes(mesh::TriMesh) = mesh.num_oshapes

import Mesh.oriented_shape_for_fe
oriented_shape_for_fe(fe::FENum, mesh::TriMesh) = mesh.fes[fe].oshapenum

import Mesh.num_side_faces_for_fe
num_side_faces_for_fe(fe::FENum, mesh::TriMesh) = ref_tri_for_fe(fe, mesh).num_side_faces

import Mesh.num_side_faces_for_shape
num_side_faces_for_shape(oshapenum::OShapeNum, mesh::TriMesh) = mesh.oshapes[oshapenum].num_side_faces

import Mesh.dependent_dim_for_oshape_side
dependent_dim_for_oshape_side(oshapenum::OShapeNum, sf::FEFaceNum, mesh::TriMesh) =
  mesh.oshapes[oshapenum].dep_dims_by_sideface[sf]

import Mesh.fe_inclusions_of_nb_side
fe_inclusions_of_nb_side(nbsn::NBSideNum, mesh::TriMesh) = mesh.nbsideincls_by_nbsidenum[nbsn]

import Mesh.nb_side_num_for_fe_side
nb_side_num_for_fe_side(fe::FENum, sf::FEFaceNum, mesh::TriMesh) = mesh.nbsidenums_by_feface[(fe, sf)]

import Mesh.is_boundary_side
is_boundary_side(fe::FENum, face::FEFaceNum, mesh::TriMesh) = !haskey(mesh.nbsidenums_by_feface, (fe,face))

import Mesh.num_boundary_sides
num_boundary_sides(mesh::TriMesh) = mesh.num_b_sides

import Mesh.shape_diameter_inv
shape_diameter_inv(oshapenum::OShapeNum, mesh::TriMesh) = mesh.oshapes[oshapenum].diameter_inv

import Mesh.max_fe_diameter
max_fe_diameter(mesh::TriMesh) = mapreduce(rt -> 1./rt.diameter_inv, max, mesh.oshapes)

import Mesh.fe_interior_origin!
function fe_interior_origin!(fe::FENum, fill::Vector{R}, mesh::TriMesh)
  const el_tri = mesh.fes[fe]
  fill[1] = el_tri.v1._1
  fill[2] = el_tri.v1._2
end


# Integration Functions

# Integrate a monomial on the indicated face of a finite element of given oriented
# shape, with the monomial interpreted locally on the face.
import Mesh.integral_face_rel_on_oshape_face
function integral_face_rel_on_oshape_face(mon::Monomial,
                                          oshapenum::OShapeNum, face::FEFaceNum,
                                          mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshapenum]

  if face == Mesh.interior_face
    # Order the interior-relative vertex points by their first coordinate values.
    const pts = sort!([zeroPt, ref_tri.v12, ref_tri.v13])
    const slope_12, slope_13 = slope_between(pts[1], pts[2]), slope_between(pts[1], pts[3])

    if pts[1]._1 == pts[2]._1 # vertical side on left between points 1 and 2
      # Integrate over the area bounded above and below by the lines between points 1 and 3 and 2 and 3,
      # and horizontally between the vertical left side formed by points 1 and 2, and point 3 on the right
      # where the upper and lower bounding lines meet.
      const slope_23 = slope_between(pts[2], pts[3])
      integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[3], slope_13, slope_23, pts[1]._1)
    else
      # Points 1 and 2 do not form a vertical line. Integrate between points 1 and 2, and between points 2 and 3 if
      # points 2 and 3 don't lie on a vertical line.
      const fst_seg = integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[1], slope_12, slope_13, pts[2]._1)
      const snd_seg =
        if pts[2]._1 == pts[3]._1 # vertical right side, no second segment
          zeroR
        else
          const slope_23 = slope_between(pts[3], pts[2])
          integral_mon_between_lines_meeting_at_point_and_vert_line(mon, pts[3], slope_13, slope_23, pts[2]._1)
        end
      fst_seg + snd_seg
    end
  else # line integral along side face
    # We want to compute int_{t=0..1} mon(p(t)-o) |p'(t)| dt, where p:[0,1] -> R^2
    # is a bijection traversing the side face smoothly, p(0) and p(1) being the side
    # endpoints, and o is the local origin for the side, which is the side's midpoint.
    const a,b = sideface_endpoint_pair(face,
                                       zeroPt, ref_tri.v12, ref_tri.v13,
                                       ref_tri.nums_faces_between_vertexes,
                                       false) # lesser endpoints first => false
    # The local origin of each side face is the midpoint.
    const o = 1/2*(a + b)
    # Path traversing the side from a to b, in local coordinates, relative to the side's local origin.
    const path_orel = let t = MON_VAR_1D
      [a._1 + (b._1 - a._1)*t - o._1,
       a._2 + (b._2 - a._2)*t - o._2]
    end
    const pullback_poly_1d = Poly.precompose_with_poly_path(mon, path_orel) * norm(b-a)
    const pullback_antider = Poly.antideriv(dim(1), pullback_poly_1d)
    Poly.value_at(pullback_antider, oneR) - Poly.value_at(pullback_antider, zeroR)
  end
end


# Integrate a global function f on the indicated finite element face, multiplied against
# a monomial interpreted locally on the face.
import Mesh.integral_global_x_face_rel_on_fe_face
function integral_global_x_face_rel_on_fe_face(f::Function,
                                               mon::Monomial,
                                               fe::FENum, face::FEFaceNum,
                                               mesh::TriMesh)
  const ref_tri = ref_tri_for_fe(fe, mesh)
  const fe_v1 = mesh.fes[fe].v1

  if face == Mesh.interior_face

    function f_x_mon(int_rel_x::Vector{R})
      const irel_1, irel_2 = int_rel_x[1], int_rel_x[2]
      const mon_val = Poly.value_at(mon, int_rel_x)
      try
        const x = int_rel_x # temporarily reuse storage of int_rel_x for global x
        # vertex 1 is the local origin for the interior
        x[1] += fe_v1._1
        x[2] += fe_v1._2
        f(x) * mon_val
      finally
        # restore int_rel_x to original contents
        int_rel_x[1] = irel_1; int_rel_x[2] = irel_2
      end
    end

    # Order the interior-relative vertex points by their first coordinate values.
    const pts = sort!([zeroPt, ref_tri.v12, ref_tri.v13])
    const slope_12, slope_13 = slope_between(pts[1], pts[2]), slope_between(pts[1], pts[3])

    if pts[1]._1 == pts[2]._1 # vertical side on left between points 1 and 2
      # Integrate over the area bounded above and below by the lines between points 1 and 3 and 2 and 3,
      # and horizontally between the vertical left side formed by points 1 and 2, and point 3 on the right
      # where the upper and lower bounding lines meet.
      const slope_23 = slope_between(pts[2], pts[3])
      integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[3], slope_13, slope_23, pts[1]._1, mesh)
    else
      # Points 1 and 2 do not form a vertical line. Integrate between points 1 and 2, and between points 2 and 3 if
      # points 2 and 3 don't lie on a vertical line.
      const fst_seg = integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[1], slope_12, slope_13, pts[2]._1, mesh)
      const snd_seg =
        if pts[2]._1 == pts[3]._1 # vertical right side, no second segment
          zeroR
        else
          const slope_23 = slope_between(pts[3], pts[2])
          integral_fn_between_lines_meeting_at_point_and_vert_line(f_x_mon, pts[3], slope_13, slope_23, pts[2]._1, mesh)
        end
      fst_seg + snd_seg
    end
  else # side face
    # We want to compute int_{t=0..1} f(p(t)) mon(p(t)-o) |p'(t)| dt, where p:[0,1] -> R^2
    # is a bijection traversing the side face smoothly, p(0) and p(1) being the side
    # endpoints, and o is the local origin for the side, which is the side's midpoint.
    const a,b = sideface_endpoint_pair(face,
                                       fe_v1, fe_v1 + ref_tri.v12, fe_v1 + ref_tri.v13,
                                       ref_tri.nums_faces_between_vertexes,
                                       false) # lesser endpoints first => false
    const o = 1/2*(a + b) # The local origin of each side face is the midpoint.
    const side_len = hypot(b._1 - a._1, b._2 - a._2)

    const integrand = let p_t = mesh.integrand_work_array
      function(t::Vector{R}) # compute f(p(t)) mon(p(t)-o) |p'(t)|
        const t = t[1]
        p_t[1] = a._1 + (b._1 - a._1)*t
        p_t[2] = a._2 + (b._2 - a._2)*t
        const f_val = f(p_t)
        # Make p_t side origin relative for the monomial evaluation.
        p_t[1] -= o._1; p_t[2] -= o._2
        const mon_val = Poly.value_at(mon, p_t)
        f_val * mon_val * side_len  # |p'(t)| = |b-a| = side length
      end
    end

    hcubature(integrand,
              singleZero, singleOne,
              mesh.integration_rel_err, mesh.integration_abs_err)[1]
  end
end

# Integrate a side-local monomial vs. a vector monomial interpreted relative to the
# entire finite element, dot multiplied with the outward normal for the side.
import Mesh.integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side
function integral_side_rel_x_fe_rel_vs_outward_normal_on_oshape_side(mon::Monomial,
                                                                     vmon::VectorMonomial,
                                                                     oshapenum::OShapeNum, side_face::FEFaceNum,
                                                                     mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshapenum]

  # fe-relative side endpoints
  const a,b = sideface_endpoint_pair(side_face,
                                     zeroPt, ref_tri.v12, ref_tri.v13,
                                     ref_tri.nums_faces_between_vertexes,
                                     false) # lesser endpoints first => false
  # The local origin of each side face is the midpoint.
  const side_o = 1/2*(a + b)
  const onormal = ref_tri.outward_normals_by_sideface[side_face]
  const side_len = hypot(b._1 - a._1, b._2 - a._2)

  # For path p:[0,1] -> S traversing from a to b linearly in fe-relative coordinates, compute
  #   mon(p(t)-side_o) vmon(p(t)).n  |p'(t)|
  const integrand = let p_t = mesh.integrand_work_array
    function(t::Vector{R})
      const t = t[1]
      # fe-relative path point
      p_t[1] = a._1 + (b._1 - a._1)*t
      p_t[2] = a._2 + (b._2 - a._2)*t
      const vmon_dot_normal = let vmon_mon_val = Poly.value_at(vmon.mon, p_t)
        vmon.mon_pos == 1 ? onormal._1 * vmon_mon_val: onormal._2 * vmon_mon_val
      end
      # Make p_t relative to the side origin for the monomial evaluation.
      p_t[1] -= side_o._1; p_t[2] -= side_o._2
      const mon_val = Poly.value_at(mon, p_t)
      mon_val * vmon_dot_normal * side_len
    end
  end

  hcubature(integrand,
            singleZero, singleOne,
            mesh.integration_rel_err, mesh.integration_abs_err)[1]
end

# Integrate a finite element relative monomial vs. a side relative monomial.
import Mesh.integral_fe_rel_x_side_rel_on_oshape_side
function integral_fe_rel_x_side_rel_on_oshape_side(fe_rel_mon::Monomial,
                                                   side_rel_mon::Monomial,
                                                   oshapenum::OShapeNum, side_face::FEFaceNum,
                                                   mesh::TriMesh)
  const ref_tri = mesh.oshapes[oshapenum]

  # fe-relative side endpoints
  const a,b = sideface_endpoint_pair(side_face,
                                     zeroPt, ref_tri.v12, ref_tri.v13,
                                     ref_tri.nums_faces_between_vertexes,
                                     false) # lesser endpoints first => false
  # The local origin of each side face is the midpoint.
  const side_o = 1/2*(a + b)
  const onormal = ref_tri.outward_normals_by_sideface[side_face]
  const side_len = hypot(b._1 - a._1, b._2 - a._2)

  # For path p:[0,1] -> S traversing from a to b linearly in fe-relative coordinates, compute
  #   side_rel_mon(p(t)-side_o) fe_rel_mon(p(t))  |p'(t)|
  const integrand = let p_t = mesh.integrand_work_array
    function(t::Vector{R})
      const t = t[1]
      # fe-relative path point
      p_t[1] = a._1 + (b._1 - a._1)*t
      p_t[2] = a._2 + (b._2 - a._2)*t
      const fe_rel_mon_val = Poly.value_at(fe_rel_mon, p_t)
      # Make p_t relative to the side origin for the side-local monomial evaluation.
      p_t[1] -= side_o._1; p_t[2] -= side_o._2
      const side_rel_mon_val = Poly.value_at(side_rel_mon, p_t)
      fe_rel_mon_val * side_rel_mon_val * side_len
    end
  end

  hcubature(integrand,
            singleZero, singleOne,
            mesh.integration_rel_err, mesh.integration_abs_err)[1]
end

# required AbstractMesh functions
####################################################################

# Integrate a monomial over the triangular region bounded on two sides by two
# non-vertical lines of indicated slopes which meet at a point p, and by the
# indicated vertical line as the remaining side.
function integral_mon_between_lines_meeting_at_point_and_vert_line(mon::Monomial,
                                                                   p::Point,
                                                                   slope_1::R, slope_2::R,
                                                                   vert_line_x::R)
  # slopes of the lower and upper lines
  const m_l, m_u = p._1 < vert_line_x ? (min(slope_1, slope_2), max(slope_1, slope_2)) :
                                        (max(slope_1, slope_2), min(slope_1, slope_2))
  # polynomials representing the lower and upper lines for the inner (y) integral's bounds
  const y_l, y_u = let x = MON_VAR_1D
    p._2 + m_l*(x-p._1),
    p._2 + m_u*(x-p._1)
  end

  # Compute inner integral, which is an integral over y with x held constant, as a polynomial in x.
  const mon_antider = Poly.antideriv(dim(2), mon)
  const inner_int = Poly.reduce_dim_by_subst(dim(2), y_u, mon_antider) - Poly.reduce_dim_by_subst(dim(2), y_l, mon_antider)

  # Difference of anti-derivative of inner integral over 1st coordinate bounds is the final integral value.
  const antider_inner_int = Poly.antideriv(dim(1), inner_int)
  # bounds of integration for outer integral
  const ib_l, ib_u = p._1 < vert_line_x ? (p._1, vert_line_x) : (vert_line_x, p._1)
  Poly.value_at(antider_inner_int, ib_u) - Poly.value_at(antider_inner_int, ib_l)
end

# Integrate a function over the triangular region bounded on two sides by two
# non-vertical lines of indicated slopes which meet at a point p, and by the
# indicated vertical line as the remaining side.
function integral_fn_between_lines_meeting_at_point_and_vert_line(g::Function,
                                                                  p::Point,
                                                                  slope_1::R, slope_2::R,
                                                                  vert_line_x::R,
                                                                  mesh::TriMesh)
  const xmin, xmax = min(p._1, vert_line_x), max(p._1, vert_line_x)
  const w = xmax - xmin # width of triangular integration region

  # Method
  # We will pull back the integration over the original triangular section T to an integration
  # over the unit square, by change of variables via the bijection
  #   t:[0,1]^2 -> T
  #   t(x,y) = (xmin + x w,  p_2 + m1 (xmin + x w - p_1) + y(m2 - m1)(xmin + x w - p_1))
  #     for (x,y) in [0,1]^2.  Here m1 and m2 are the slopes of the non-vertical bounding lines.
  # The determinant of the derivative matrix Dt(x,y) is
  #                     |          w                           0               |
  #  det Dt(x,y)| = det |                                                      |
  #                     | m1 w + y(m2 - m1) w      (m2 - m1)(xmin + x w - p_1) |
  #               = w (m2 - m1) (xmin + x w - p_1).
  # Now by applying change of variables in the integral of g over the T via the mapping t,
  # we have
  #   int_T g = int_0^1 int_0^1 g(t(x,y)) |det Dt(x,y)| dy dx
  #            = int_0^1 int_0^1 g(t(x,y)) |w (m2 - m1) (xmin + x w - p_1)| dy dx

  const slopediff = slope_2 - slope_1
  const w_slopediff = w * slopediff

  function gt_absdetDt(r::Vector{R}) # evaluate integrand at rectangle cubature point r
    const x, y = r[1], r[2]
    const xT = xmin + x * w      # (triangle x)
    const xT_minus_p1 = xT - p._1
    try
      # Construct the triangle point t(x,y).
      const t = r # reuse the rectangle point array temporarily for the triangle point
      t[1] = xT
      t[2] = p._2 + slope_1 * xT_minus_p1  +  y * slopediff * xT_minus_p1
      const det_Dt = w_slopediff * xT_minus_p1
      g(t) * abs(det_Dt)
    finally
      # restore input array to its original state
      r[1] = x; r[2] = y
    end
  end

  hcubature(gt_absdetDt,
            mesh.space_dim_zeros, mesh.space_dim_ones, # unit square bounds
            mesh.integration_rel_err, mesh.integration_abs_err)[1]
end


####################################################################
# Mesh Construction

# Construct from Gmsh .msh formatted input stream, and number of subdivision operations to perform within each
# mesh element from the stream.  Some elements may have additional iterations specified via Gmsh tags.
function TriMesh(ios::IO, base_subdiv_iters::Integer, integration_rel_err::R = 1e-12, integration_abs_err::R = 1e-12)

  const base_subdiv_iters = base_subdiv_iters >= 0 ? uint(base_subdiv_iters) : throw(ArgumentError("Number of iterations must be non-negative."))

  # Read the points from the input mesh, storing by mesh node number.  The elements section of the input
  # stream will refer to these node numbers.
  const mesh_pts_by_nodenum = read_points(ios)

  # Read to the beginning of the elements section.
  if !read_through_line(ios, "\$Elements\n") error("Could not find beginning of elements section in mesh file.") end

  # Skip point and line elements, estimating number of input mesh triangles from the number of elements declared.
  const decl_num_els = uint64(readline(ios)) # This count can include unwanted lower order elements.
  line, lines_read = read_until(ios, is_polytope_or_endmarker)
  const est_mesh_tris = decl_num_els - (lines_read - 1)

  # Make estimates of the number of finite elements and reference triangles for storage allocation.
  # This estimate ignores optional additional subdivisions that may be specified for some input elements.
  const est_fes = est_mesh_tris * 4^base_subdiv_iters
  # We'll estimate 2 reference triangles for each mesh element if we are subdividing, 1 if not. We ignore
  # for this estimate any additional reference triangles for hanging node support.
  const est_ref_tris = (base_subdiv_iters > 0 ? 2 : 1) * est_mesh_tris
  println("Estimating $(int(est_fes)) finite elements, $(int(est_ref_tris)) reference triangles.")

  # Data to be updated as input element lines are processed.
  const oshapes = sizehint(Array(RefTri,0), est_ref_tris)
  const fes = sizehint(Array(ElTri, 0), est_fes)
  const fe_faces_by_endpoints = sizehint(Dict{(Point,Point), Array{(FENum,FEFaceNum),1}}(), uint64(est_fes * 3/2)) # estimate assumes most fes have 3 sides, each in 2 fes
  num_nb_sides = 0
  num_b_sides = 0

  # Function via which all finite elements will be created during the processing of input element lines.
  const fe_maker = function(oshapenum::OShapeNum, v1::Point, v2::Point, v3::Point)
     push!(fes, ElTri(oshapenum, v1))
     const fe_num = fenum(length(fes))
     const nums_faces_btw_verts = oshapes[oshapenum].nums_faces_between_vertexes
     const nb_delta, b_delta = register_fe_faces_by_endpoints(fe_num, v1,v2,v3, nums_faces_btw_verts, fe_faces_by_endpoints)
     num_nb_sides += nb_delta
     num_b_sides += b_delta
  end

  # process input mesh elements
  while line != "\$EndElements\n" && line != ""
    process_el_line(line, mesh_pts_by_nodenum, base_subdiv_iters, fe_maker, oshapes)
    line = readline(ios)
  end

  if line != "\$EndElements\n" error("Missing end of elements marker in mesh file.") end

  # Create the final non-boundary sides data structures based on the mapping of side endpoints to fe faces.
  const nbsidenums_by_feface, nbsideincls_by_nbsidenum  = create_nb_sides_data(fe_faces_by_endpoints, num_nb_sides_hint=num_nb_sides)

  const space_dim = dim(2)

  TriMesh(fes,
          oshapes,
          nbsidenums_by_feface,
          nbsideincls_by_nbsidenum,
          num_b_sides,
          space_dim,
          Monomial(zeros(Deg,2)),
          fenum(length(fes)),
          nbsidenum(length(nbsideincls_by_nbsidenum)),
          Mesh.oshapenum(length(oshapes)),
          Array(R, space_dim), # integrand_work_array
          zeros(R, space_dim),
          ones(R,  space_dim),
          integration_rel_err,
          integration_abs_err)

end


function process_el_line(line::String,
                         mesh_pts_by_nodenum::Array{Point,1},
                         base_subdiv_iters::Uint,
                         fe_maker::Function,
                         oshapes::Array{RefTri,1}) # will be appended to
  const toks = split(line, ' ')

  const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
  if is_lower_order_el_type(el_type) return end
  if !is_3_node_triangle_el_type(el_type) error("Element type $el_type is not supported.") end

  const el_tags = toks[TOKN_ELLINE_NUMTAGS+1 : TOKN_ELLINE_NUMTAGS+int(toks[TOKN_ELLINE_NUMTAGS])]

  const v1,v2,v3 = lookup_vertexes(toks, mesh_pts_by_nodenum)

  # Add any extra subdivision iterations to be done in this element.
  const subdiv_iters = base_subdiv_iters + extra_subdiv_iters(v1,v2,v3, el_tags)

  # For each of the three sides of this mesh element to be subdivided, the generated subdivision
  # elements can be made to have more than one face between vertexes lying on the side.  This is
  # used to support "hanging" nodes where a finer subdivision is adjacent to this one.  The triplet
  # returned represents the number of faces between element vertex pairs embedded in the input mesh
  # element sides from v1 to v2, v2 to v3, and v3 to v1, respectively.
  const nums_faces_btw_verts = nums_faces_between_vertexes(v1,v2,v3, el_tags)

  # Register the primary reference triangles for our mesh element's subdivisions.
  const pri_oshapenums_by_nums_faces_btw_verts = register_primary_ref_tris(v1,v2,v3, nums_faces_btw_verts, subdiv_iters, oshapes)

  if subdiv_iters > 0
    # Register the secondary reference triangle.
    push!(oshapes, secondary_ref_tri(v1,v2,v3, subdiv_iters))
    const sec_oshapenum = Mesh.oshapenum(length(oshapes))

    # Do the subdivisions.
    subdivide_primary(v1,v2,v3,
                      subdiv_iters,
                      nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)
  else # no subdivision to be done
    # The mesh element itself is our finite element, with its own reference triangle (added in register_primary_ref_tris above).
    const oshapenum = pri_oshapenums_by_nums_faces_btw_verts(nums_faces_btw_verts)
    fe_maker(oshapenum, v1,v2,v3)
  end
end # process_el_line


# Create primary reference triangles for the given triangle to be subdivided. These reference triangles
# are rescaled translations of the original undivided triangle. If nums_faces_btw_verts is other than
# (1,1,1), then multiple reference triangles will be generated, differing in the number of side faces
# between their vertexes for support of "hanging" nodes.  The function returns a function mapping the
# numbers of side faces between vertexes to the newly registered reference triangle (oshape) number.
# Pains are taken to not create more reference triangles than are actually used, so that some global
# mesh properties such as maximum element diameter can be determined from only the reference elements.
function register_primary_ref_tris(v1::Point, v2::Point, v3::Point,
                                   nums_faces_btw_verts::(Int,Int,Int),
                                   subdiv_iters::Integer,
                                   oshapes::Array{RefTri,1})
  const orig_num_oshapes = length(oshapes)

  if subdiv_iters == 0
    # Not subdividing, just add the passed triangle itself as the only reference triangle.
    push!(oshapes, primary_ref_tri(v1,v2,v3, nums_faces_btw_verts, 0))

  else # at least one subdivision
    # We need to add a separate primary reference triangle for each triplet of number of faces between
    # vertexes that may occur in primary triangles in the subdivision. We use the set of triplets below
    # to record the various numbers of faces per side required for the reference triangles prior to
    # creating them.
    const nums_sides_btw_verts_triplets = Set{(Int,Int,Int)}()

    # To support the corner elements, we need primary reference triangles which get 2 of their 3 numbers
    # of faces between vertex pairs from nums_faces_btw_verts, because they have two sides on the original
    # undivided triangle's sides. Note that if we're just doing one subdivision then each of the three
    # primary triangles produced are corner triangles so we're done.
    for i=1:3
      const j = mod1(i+1,3) # next vertex pair number, cyclically
      add!(nums_sides_btw_verts_triplets,
           (i==1||j==1 ? nums_faces_btw_verts[1] : 1,
            i==2||j==2 ? nums_faces_btw_verts[2] : 1,
            i==3||j==3 ? nums_faces_btw_verts[3] : 1))
    end

    # If subdividing more than once, then we also need:
    #   1) a primary reference triangle with only one face between its vertex pairs (always present within
    #      the secondary triangle produced in the first subdivision), and
    #   2) primary reference triangles for each choice of one number of faces from nums_faces_btw_verts for
    #      some pair of vertexes, and just one face between the other two vertex pairs.
    if subdiv_iters >= 2
      add!(nums_sides_btw_verts_triplets, ONE_FACE_BETWEEN_VERTEX_PAIRS)
      for i=1:3
        add!(nums_sides_btw_verts_triplets,
             (i==1 ? nums_faces_btw_verts[1] : 1,
              i==2 ? nums_faces_btw_verts[2] : 1,
              i==3 ? nums_faces_btw_verts[3] : 1))
      end
    end

    for faces_per_vert_pair in nums_sides_btw_verts_triplets
      push!(oshapes, primary_ref_tri(v1,v2,v3, faces_per_vert_pair, subdiv_iters))
    end
  end

  # Return a lookup function for these new primary oshape numbers by numbers of side faces between vertexes.
  const first_new_oshapenum = Mesh.oshapenum(orig_num_oshapes + 1)
  const last_new_oshapenum = Mesh.oshapenum(length(oshapes))
  function(nums_faces_btw_verts::(Int,Int,Int))
    for osn=first_new_oshapenum:last_new_oshapenum
      if oshapes[osn].nums_faces_between_vertexes == nums_faces_btw_verts
        return osn
      end
    end
    error("Reference triangle not found by numbers of side faces between vertexes: $nums_faces_btw_verts.")
  end
end

function primary_ref_tri(v1::Point, v2::Point, v3::Point,
                         nums_faces_btw_verts::(Int,Int,Int),
                         subdiv_iters::Integer)
  const scale = 1./2^subdiv_iters
  RefTri(v1,v2,v3, scale, nums_faces_btw_verts)
end

# Returns the reference element for the secondary, less numerous, inverted triangles which
# form in the middle positions after triangle subdivisions. Note that passed vertexes represent
# those of an outermost triangle to be subdivided, and hence describe a triangle of the primary,
# not secondary, oriented shape.
function secondary_ref_tri(pri_v1::Point, pri_v2::Point, pri_v3::Point, # vertexes of the undivided, *primary* shape
                           subdiv_iters::Integer)
  if subdiv_iters == 0 error("Secondary reference triangle requires non-zero number of iterations.") end
  const v1,v2,v3 = .5(pri_v1+pri_v2), .5(pri_v2+pri_v3), .5(pri_v3+pri_v1)
  const scale = 1./2^(subdiv_iters-1) # v1,v2,v3 already represents one subdivision
  RefTri(v1,v2,v3, scale, ONE_FACE_BETWEEN_VERTEX_PAIRS)
end


function subdivide_primary(v1::Point, v2::Point, v3::Point,
                           iters::Uint,
                           nums_faces_btw_verts::(Int,Int,Int),
                           pri_oshapenums_by_nums_faces_btw_verts::Function, sec_oshapenum::OShapeNum,
                           fe_maker::Function)
  if iters == 0
    const oshapenum = pri_oshapenums_by_nums_faces_btw_verts(nums_faces_btw_verts)
    fe_maker(oshapenum, v1,v2,v3)
  else
    const midpt_12 = 0.5(v1+v2)
    const midpt_13 = 0.5(v1+v3)
    const midpt_23 = 0.5(v2+v3)

    # The sub-triangle including v1 has its first and last vertex pairs embedded in the original triangle's first and
    # third sides, and so inherit the corresponding numbers of faces between vertexes from nums_faces_btw_verts.
    # The middle vertex pair (midpt_12 to midpt_13) of this sub-triangle does not lie along an original side, so has
    # only one face between these vertexes.  Similar arguments apply for the remaining sub-triangles.
    const st1_nums_faces_btw_verts = (nums_faces_btw_verts[1], 1, nums_faces_btw_verts[3])
    subdivide_primary(v1, midpt_12, midpt_13,
                      iters-1,
                      st1_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    const st2_nums_faces_btw_verts = (nums_faces_btw_verts[1], nums_faces_btw_verts[2], 1)
    subdivide_primary(midpt_12, v2, midpt_23,
                      iters-1,
                      st2_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    const st3_nums_faces_btw_verts = (1, nums_faces_btw_verts[2], nums_faces_btw_verts[3])
    subdivide_primary(midpt_13, midpt_23, v3,
                      iters-1,
                      st3_nums_faces_btw_verts,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)

    # secondary sub-triangle
    subdivide_secondary(midpt_12, midpt_23, midpt_13,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)
  end
end


function subdivide_secondary(v1::Point, v2::Point, v3::Point,
                             iters::Uint,
                             pri_oshapenums_by_nums_faces_btw_verts::Function, sec_oshapenum::OShapeNum,
                             fe_maker::Function)
  if iters == 0
    fe_maker(sec_oshapenum, v1,v2,v3)
  else
    const midpt_12 = 0.5(v1+v2)
    const midpt_13 = 0.5(v1+v3)
    const midpt_23 = 0.5(v2+v3)

    subdivide_secondary(v1, midpt_12, midpt_13,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    subdivide_secondary(midpt_12, v2, midpt_23,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    subdivide_secondary(midpt_13, midpt_23, v3,
                        iters-1,
                        pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                        fe_maker)

    # secondary sub-triangle
    subdivide_primary(midpt_13, midpt_12, midpt_23,
                      iters-1,
                      ONE_FACE_BETWEEN_VERTEX_PAIRS,
                      pri_oshapenums_by_nums_faces_btw_verts, sec_oshapenum,
                      fe_maker)
  end
end


function lookup_vertexes(toks::Array{String,1}, mesh_pts_by_nodenum::Array{Point,1})
  const last_tag_tokn = TOKN_ELLINE_NUMTAGS + int(toks[TOKN_ELLINE_NUMTAGS])
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 1])],
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 2])],
  mesh_pts_by_nodenum[uint64(toks[last_tag_tokn + 3])]
end


# Register the finite element's side faces by their endpoints in the passed registry.
function register_fe_faces_by_endpoints(fe_num::FENum,
                                        v1::Point, v2::Point, v3::Point,
                                        nums_faces_btw_verts::(Int,Int,Int),
                                        fe_faces_by_endpoints::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}})
  const sf_endpt_pairs = sideface_endpoint_pairs(v1,v2,v3, nums_faces_btw_verts,
                                                 true) # lesser endpoints first
  nb_sides_delta = 0
  b_sides_delta = 0

  # Record the side faces included by this elementsment.
  for sf=fefacenum(1):fefacenum(length(sf_endpt_pairs))
    const sf_endpts = sf_endpt_pairs[sf]
    const existing_side_incls = get(fe_faces_by_endpoints, sf_endpts, nothing)
    if existing_side_incls != nothing
      push!(existing_side_incls, (fe_num, sf))
      b_sides_delta -= 1
      nb_sides_delta += 1
    else
      fe_faces_by_endpoints[sf_endpts] = [(fe_num, sf)]
      b_sides_delta += 1
    end
  end
  nb_sides_delta, b_sides_delta
end

# Return the triplet consisting of
#   1) nb side numbers by (fe,face) :: Dict{(FENum, FEFaceNum), NBSideNum},
#   2) an array of NBSideInclusions indexed by nb side number.
#   3) an array of side dependent dimension by nb side number
function create_nb_sides_data(fe_faces_by_endpoints::Dict{(Point, Point), Array{(FENum, FEFaceNum),1}};
                              num_nb_sides_hint::Integer=0)
  const nbsidenums_by_feface = sizehint(Dict{(FENum, FEFaceNum), NBSideNum}(), 2*num_nb_sides_hint)
  const nbsideincls_by_nbsidenum = sizehint(Array(NBSideInclusions, 0), num_nb_sides_hint)

  for side_endpoints in keys(fe_faces_by_endpoints)
    const fe_faces = fe_faces_by_endpoints[side_endpoints]
    const num_fe_incls = length(fe_faces)
    if num_fe_incls == 2
      sort!(fe_faces)
      const (fe_1, sf_1) = fe_faces[1]
      const (fe_2, sf_2) = fe_faces[2]

      # Assign a new non-boundary side number.
      const nb_side_num = nbsidenum(length(nbsideincls_by_nbsidenum) + 1)

      # Add nb side inclusions structure for this nb side number.
      push!(nbsideincls_by_nbsidenum,
            NBSideInclusions(nb_side_num,
                             fe_1, sf_1,
                             fe_2, sf_2))

      # Map the two fe/face pairs to this nb side number.
      nbsidenums_by_feface[fe_faces[1]] = nb_side_num
      nbsidenums_by_feface[fe_faces[2]] = nb_side_num
    elseif num_fe_incls > 2
      error("Invalid mesh: side encountered with more than two including finite elements.")
    end
  end

  nbsidenums_by_feface, nbsideincls_by_nbsidenum
end

# Computing Outward Normals for Side Faces
# The outward normal for an inter-vertex vector w (extended to R^3 via 0 third component) is
#   n = (w x (0,0,cc)) / |w|
#     = cc (w[2], -w[1], 0) / |w|.
# where cc = 1 if w is part of a clockwise traversal of the triangle, and -1 otherwise.
# For any pair of successive inter-vertex vectors w_1 and w_2, cc can be computed as:
#   cc = sgn(z(w_1 x w_2)), where z is projection of component 3.
function outward_normals(v1::Point, v2::Point, v3::Point, nums_faces_btw_verts::(Int,Int,Int)=(1,1,1))
  const normals = sizehint(Array(Vec, 0), sum(nums_faces_btw_verts))
  # inter-vertex vectors
  const v12 = v2-v1
  const v23 = v3-v2
  const v31 = v1-v3
  const cc = counterclockwise(v12, v23) ? 1. : -1.
  for (ivv, num_faces_btw_verts) in ((v12, nums_faces_btw_verts[1]),
                                     (v23, nums_faces_btw_verts[2]),
                                     (v31, nums_faces_btw_verts[3]))
    const len = norm(ivv)
    const n = Vec(cc * ivv._2/len, -cc * ivv._1/len)
    for sf=1:num_faces_btw_verts
      push!(normals, n)
    end
  end
  normals
end

function sideface_endpoint_pair(sf::FEFaceNum,
                                v1::Point, v2::Point, v3::Point,
                                nums_faces_btw_verts::(Int,Int,Int),
                                lesser_endpts_first::Bool=false)
  cur_sf = Mesh.fefacenum(0)
  const vert_pairs = ((v1,v2), (v2,v3), (v3,v1))
  for i=1:3
    const va,vb = vert_pairs[i]
    const num_faces_btw_va_vb = nums_faces_btw_verts[i]
    if num_faces_btw_va_vb == 1
      cur_sf += 1
      if cur_sf == sf
        return mk_endpoint_pair(va, vb, lesser_endpts_first)
      end
    elseif num_faces_btw_va_vb == 2
      cur_sf += 1
      if cur_sf == sf
        const midpt = 0.5(va+vb)
        return mk_endpoint_pair(va, midpt, lesser_endpts_first)
      end
      cur_sf += 1
      if cur_sf == sf
        const midpt = 0.5(va+vb)
        return mk_endpoint_pair(midpt, vb, lesser_endpts_first)
      end
    else
      error("Only 1 or 2 faces between triangle vertexes are currently supported.")
    end
  end
  error("Side face $sf not found.")
end

# This function defines the side faces enumeration for a finite element of given vertexes and numbers
# of faces between vertexes. Side faces are returned as an array of side endpoint pairs indexed by
# side face number. If lesser_endpts_first is true, then each endpoint pair endpoint will have
# the lesser point (compared lexicographically) in the first component of the pair.
function sideface_endpoint_pairs(v1::Point, v2::Point, v3::Point,
                                 nums_faces_btw_verts::(Int,Int,Int),
                                 lesser_endpts_first::Bool)
  const num_side_faces = sum(nums_faces_btw_verts)
  const sf_endpt_pairs = sizehint(Array((Point,Point),0), num_side_faces)
  const vert_pairs = ((v1,v2), (v2,v3), (v3,v1))
  for i=1:3
    const va,vb = vert_pairs[i]
    const num_faces_btw_va_vb = nums_faces_btw_verts[i]
    if num_faces_btw_va_vb == 1
      push!(sf_endpt_pairs, mk_endpoint_pair(va, vb, lesser_endpts_first))
    elseif num_faces_btw_va_vb == 2
      const midpt = 0.5(va+vb)
      push!(sf_endpt_pairs, mk_endpoint_pair(va, midpt, lesser_endpts_first))
      push!(sf_endpt_pairs, mk_endpoint_pair(midpt, vb, lesser_endpts_first))
    else
      error("Only 1 or 2 faces between triangle vertexes are currently supported.")
    end
  end
  sf_endpt_pairs
end

function mk_endpoint_pair(pt1::Point, pt2::Point, lesser_pt_first::Bool)
  if lesser_pt_first isless(pt1, pt2) ? (pt1, pt2) : (pt2, pt1)
  else (pt1, pt2) end
end


function sideface_dep_dims(v1::Point, v2::Point, v3::Point, nums_faces_btw_verts::(Int,Int,Int))
  const dep_dims = sizehint(Array(Dim,0), sum(nums_faces_btw_verts))
  const vert_pairs = ((v1,v2), (v2,v3), (v3,v1))
  for i=1:3
    const va,vb = vert_pairs[i]
    const ddim = side_dep_dim(va, vb)
    for i=1:nums_faces_btw_verts[i]
      push!(dep_dims, ddim)
    end
  end
  dep_dims
end

# Return a dependent dimension for the given side endpoints.
side_dep_dim(ep1::Point, ep2::Point) = abs(ep2._1-ep1._1) >= abs(ep2._2-ep1._2) ? dim(2) : dim(1)


function extra_subdiv_iters(v1::Point, v2::Point, v3::Point, el_tags::Array{String,1})
  # TODO: Allow specifying extra iterations by tagging a mesh element.
  0
end

# For each of the three sides of the indicated mesh element to be subdivided, we can specify here the
# number of faces that the generated subtriangles should have between their vertexes which lie on this side.
# This allows these subdivision elements to meet those of a finer subdivision in a mesh element adjacent
# to this element ("hanging nodes").  The numbers are returned as a triplet of integers, corresponding to
# the face counts to be generated for elements' sides within sides v1v2, v2v3, and v3v1.
function nums_faces_between_vertexes(v1::Point, v2::Point, v3::Point, el_tags::Array{String,1})
  # TODO: Allow specifying somewhere such as in mesh element tags.
  ONE_FACE_BETWEEN_VERTEX_PAIRS
end

# file reading utilities

function read_points(ios::IO)
  if !read_through_line(ios, "\$Nodes\n")
    error("Nodes section not found in mesh input file.")
  else
    # Next line should be the vert count.
    const count = int(readline(ios))
    const pts = Array(Point, count)
    l = strip(readline(ios))
    while l != "\$EndNodes" && l != ""
      const toks = split(l, ' ')
      const pt = Point(convert(R, float64(toks[TOKN_NODELINE_POINT1])), convert(R, float64(toks[TOKN_NODELINE_POINT2])))
      if length(toks) >= TOKN_NODELINE_POINT3 && toks[TOKN_NODELINE_POINT3] != "0" && toks[TOKN_NODELINE_POINT3] != "0.0"
        error("Nodes with non-zero third coordinates are not supported in this 2D mesh reader.")
      end
      pts[uint64(toks[TOKN_NODELINE_ELNUM])] = pt
      l = strip(readline(ios))
    end
    if l == "" error("End of nodes section not found in mesh file.") end
    pts
  end
end

function is_polytope_or_endmarker(line::ASCIIString)
  if line == "\$EndElements"
    true
  else
    const toks = split(line, ' ')
    !is_lower_order_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
  end
end

function read_through_line(io::IO, line::ASCIIString)
  l = readline(io)
  while l != "" && l != line
    l = readline(io)
  end
  l != ""
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IO, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end



# definitions related to the Gmsh format

typealias ElTypeNum Uint8
eltypenum(i::Integer) = if i>=1 && i<=typemax(ElTypeNum) convert(ElTypeNum, i) else error("invalid element type: $i") end
eltypenum(s::String) = eltypenum(int(s))

# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)

function is_3_node_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3

# Mesh Construction
####################################################################


###############
# Etc

const ONE_FACE_BETWEEN_VERTEX_PAIRS = (1,1,1)

const MON_VAR_1D = Monomial(1)

# Reference triangle for a finite element.
ref_tri_for_fe(fe::FENum, mesh::TriMesh) = mesh.oshapes[mesh.fes[fe].oshapenum]

# vector operations

import Base.(+)
+(v::Vec, w::Vec) = Vec(v._1 + w._1, v._2 + w._2)

import Base.(-)
-(v::Vec, w::Vec) = Vec(v._1 - w._1, v._2 - w._2)

import Base.(*)
*(r::Real, v::Vec) = Vec(r*v._1, r*v._2)

import Base.norm
norm(v::Vec) = hypot(v._1, v._2)

import Base.isless
isless(v::Vec, w::Vec) = v._1 < w._1 || v._1 == w._1 && v._2 < w._2

import Base.dot
dot(v::Vec, w::Vec) = v._1 * w._1 + v._2 * w._2


counterclockwise(u::Vec, v::Vec) =
  u._1*v._2 - u._2*v._1 > 0

slope_between(p::Point, q::Point) = (q._2 - p._2)/(q._1 - p._1)

# string representations
import Base.string, Base.show, Base.print

function string(m::TriMesh)
  "TriMesh: $(length(m.fes)) elements, $(length(m.oshapes)) reference elements, $(length(m.nbsideincls_by_nbsidenum)) non-boundary sides, $(m.num_b_sides) boundary sides"
end
print(io::IO, m::TriMesh) = print(io, string(m))
show(io::IO, m::TriMesh) = print(io, m)

end # end of module
