#function gethostname()
#  hostname = Array(Uint8, 128)
#  ccall( (:gethostname, "libc"), Int32,
#        (Ptr{Uint8}, Uint),
#        hostname, length(hostname))
#  return bytestring(convert(Ptr{Uint8}, hostname))
#end


#function compute_dot(DX::Vector, DY::Vector)
#  assert(length(DX) == length(DY))
#  n = length(DX)
#  incx = incy = 1
#  product = ccall( (:ddot_, "libLAPACK"),
#                  Float64,
#                  (Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
#                  &n, DX, &incx, DY, &incy)
#  return product
#end

#const cubature_lib = dlopen("cubature")


# See https://github.com/JuliaLang/julia/issues/1096
function integrate_R2toR(f::Function, bound_mins::Vector{Float64}, bound_maxs::Vector{Float64},
                         max_evals::Uint32, req_abs_error::Float64, req_rel_error::Float64)
  local cf = cfunction(f, Float64, (Float64,Float64)),
        res = Array(Float64, 1),
        err = Array(Float64, 1)

  status = ccall((:intR2toR,"cubature"), Int32, (Ptr{None}, Ptr{Float64}, Ptr{Float64}, Uint32, Float64, Float64, Ptr{Float64}, Ptr{Float64}),
                 cf, bound_mins, bound_maxs,
                 max_evals, req_abs_error, req_rel_error,
                 res, err)
  if status == 0
    res[1]
  else
    error("integration failed with error code $status.")
  end
end
