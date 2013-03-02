

function basis_ips(basis::Array)
  local n = length(basis), a = zeros(R, n,n)
  for i = 1:n, j = i:n
    let ip = ip_L2_unit_square(basis[i], basis[j])
      a[i,j] = ip
      if j != i
        a[j,i] = ip
      end
    end
  end
  a
end

fun_ips_vs_basis(f::Function, basis::Array) =
  Float64[ip_L2_unit_square(f, e) for e in basis]


# TODO: return a proper polynomial type here?
linear_comb(coefs::Array{Float64}, basis::Array) =
  let n = length(basis)
    (x,y) -> sum( {coefs[i]*basis[i](x,y) for i in 1:n} )
  end

# Project a function onto the subspace spanned by the passed basis.
project(f::Function, basis::Array{Function}) =
  # <Proj(f),e_i> = <f,e_i> for all e_i in the passed subspace basis.
  # Since Proj(f) = sum a_j e_j for some a_j's, this gives us a linear system which we can solve for the a_j.
  let coefs = basis_ips(basis) \ fun_ips_vs_basis(f, basis)
    linear_comb(coefs, basis)
  end

ip_L2_unit_square(f1::Function, f2::Function) =
  dbl_int_unit_square((x::Float64, y::Float64) -> dot(f1(x,y), f2(x,y)))

dbl_int_unit_square(f::Function) = f(0.5,0.5) # TODO
  #Calculus.adaptive_simpsons_outer(x::Float64 -> Calculus.adaptive_simpsons_outer(y::Float64 -> f(x,y), 0.0, 1.0, 10e-12, 60), 0.0, 1.0, 10e-12, 60)
  #Calculus.integrate(x::Float64 -> Calculus.integrate(y::Float64 -> f(x,y), 0.0, 1.0), 0.0, 1.0)




# test functions
fv1(x::Real,y::Real) = [10*sin(pi * x) * cos(pi * x * y), x^5 - 3 * x^2 * y - x * y]

fs(x::Real,y::Real) = 10*sin(pi * x) * cos(pi * x * y)

function test_vec_proj()
  local vbasis3 = poly_vec_basis(3),
        vbasis5 = poly_vec_basis(5),
        fvproj3 = project(fv1, vbasis3),
        fvproj5 = project(fv1, vbasis5)
  vdiff3(x,y) = fv1(x,y) - fvproj3(x,y)
  vdiff5(x,y) = fv1(x,y) - fvproj5(x,y)
  println("Vector diffs, projection to space of polynomials:")
  for i in 1:3, j in 1:3
    println("degree 3: ", vdiff3(i/3.0,j/3.0))
    println("degree 5: ", vdiff5(i/3.0,j/3.0))
  end
end

m2d = Monomial(1,2)
vm1_2d = VectorMonomial(m2d, dim_ix(1), dim(2))
vm2_2d = VectorMonomial(m2d, dim_ix(2), dim(2))
x2d = Monomial(1,0)
y2d = Monomial(0,1)
p1_2d = x2d + y2d
p2_2d = 3x2d + -2.0y2d
