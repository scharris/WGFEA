module Common
export Deg, deg, check_degs, Dim, dim, R, zeroR, oneR, FunctionOrConst

typealias Deg Uint8
deg(k::Integer) = check_degs(k) && convert(Deg,k)
check_degs(ks...) = if any(k -> k < 0 || k > 255, ks) error("invalid degree") else true end

typealias Dim Uint8
dim(d::Integer) = if d < 0 || d > 255 error("invalid dimension") else convert(Dim, d) end

typealias R Float64
const zeroR = zero(R)
const oneR = one(R)

typealias FunctionOrConst Union(Function, R)

end # end of module
