## The utils.jl file contains utility functions for working with polynomials.

"""
    rand_poly(::Type{T}=Float64, X, d; coeffs = I -> randn(T), n = length(X))

Returns a random polynomial with coefficients of type `T` in variables `X` with degree `d`.

# Keyword arguments:
 - `coeffs`: A function that takes a tuple of exponents and returns a coefficient.
 - `n`: The number of variables.
"""
function rand_poly(::Type{T}, X, d; coeffs = I -> randn(T), n = length(X))
    sum(I -> coeffs(I)*_monomial(X, I), Combinatorics.multiexponents(n, d))
end
rand_poly(X, d; coeffs = I -> randn(Float64), n = length(X)) =
    rand_poly(Float64, X, d; coeffs = coeffs, n = n)

"""
    Id(::Type{T}, n) where T

Returns the identity matrix of size `n` and type `T`.
"""
function Id(::Type{T}, n)
    LA.Diagonal(ones(T, n))
end
Id(n) = Id(Float64, n)
# diagonal(v) = cat(v...; dims = (1, 2)) # Without LinearAlgebra.jl


"""
    _norm(v)

Returns the Euclidean norm of a vector `v`.
"""
function _norm(v)
    sqrt(sum(x -> x^2, v))
end

"""
    Wnorm_term(I, C)

Returns the term of the Weyl norm corresponding to the exponents `I` and coefficient `C`.
"""
function Wnorm_term(I, C)
    _norm(C)^2*inv(Combinatorics.multinomial(I...))
end

"""
    Wnorm_term(IC) = Wnorm_term(IC...)

"""
function Wnorm_term(IC)
    Wnorm_term(IC...)
end

"""
    Wnorm_poly(poly; vars, expanded = false)

Returns the Weyl norm of a polynomial `poly`.

# Keyword arguments:
 - `vars`: The variables of the polynomial.
 - `expanded`: A boolean indicating whether the polynomial should be expanded before calculating the Weyl norm.
"""
function Wnorm_poly(poly;
                    vars,
                    expanded = false,
                    )
    sum(Wnorm_term, _exposcoeffs(poly; vars = vars, expanded = expanded))
end

"""
    Wnorm(F; vars = _get_vars(F), expanded = true)

Returns the Weyl norm of a system of polynomials `F`.

# Keyword arguments:
 - `vars`: The variables of the system.
 - `expanded`: A boolean indicating whether the polynomials should be expanded before calculating the Weyl norm.
"""
function Wnorm(F; vars = _get_vars(F), expanded = true)
    sum(f -> Wnorm_poly(f; vars = vars, expanded = expanded), _get_polys(F))
end

"""
    Δm1(F)

Returns the inverse of the degree matrix of a system of polynomials `F`.
"""
function Δm1(F)
    inv(LA.Diagonal(_degrees(F)))
end


"""
    D(Jx, x, n = length(x))

Returns the matrix `Jx` multiplied by the difference between the identity matrix and the outer product of `x` and its conjugate transpose.

# Keyword arguments:
 - `n`: The size of the identity matrix.
"""
function D(Jx, x, n = length(x))
    Jx*(Id(n) - x*transpose(conj(x)))
end
# D(J, x, n = length(x)) = J(x)*(Id(n) - x*transpose(conj(x)))

const ε = 1e-6 # bound for meaning zero

"""
    K(Wnorm, Δm1)

Returns a function that calculates the condition number of a system of polynomials given its Weyl norm and inverse degree matrix.
"""
function K(Wnorm, Δm1)
    (xfJ) -> begin
        M =  Δm1*D(xfJ[3],xfJ[1])
        svals = LA.svdvals(M)
        σ = svals[findlast(>(ε), svals)]
        f2 = _norm(xfJ[2])^2
        return Wnorm/sqrt(f2 + σ)
    end
end

"""
    _monomial(X, I, n = length(X))

Returns the monomial of variables `X` with exponents `I`.

# Keyword arguments:
 - `n`: The number of variables.
"""
function _monomial(X, I, n=length(X))
    prod(i -> X[i]^I[i], 1:n)
end

"""
    _exposcoeffs(poly; vars, expanded = false)

Returns a list of tuples containing the exponents and coefficients of a polynomial `poly`.

# Keyword arguments:
 - `vars`: The variables of the polynomial.
 - `expanded`: A boolean indicating whether the polynomial should be expanded before extracting the exponents and coefficients.
"""
function _exposcoeffs end

"""
    _get_vars(F)

Returns the variables of a system of polynomials `F`.
"""
function _get_vars(F) end

"""
    _get_polys(F)

Returns the polynomials of a system `F`.
"""
function _get_polys end

"""
    _degrees(F)

Returns the degrees of a system of polynomials `F`.
"""
function _degrees end

"""
    K(F; compiled = true)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
function K end

# # monomial(X, I) = prod(X.^I)

# function monomial_base(X, d; affine::Bool = false)
#     n = length(X)
#     f(I) = (I, _monomial(X, I, n))
#     return Iterators.map(f, _monomials_exponents(n, d; affine))
# end


# ## Adapted from HC.jl
# function _monomials_exponents(n, d; affine::Bool)
#     if affine
#         ## Iter onver homogeneous monomials in _n+1_ vars
#         E = map(Combinatorics.multiexponents(n+1, d)) do e
#             ## Dehomoge at collect
#             e[begin+1:end]
#         end
#     else
#         E = map(Combinatorics.multiexponents(n, d)) do e
#             collect(e)
#         end
#     end
#     sort!(E, lt = ModelKit.td_order)
#     E
# end
# function ModelKit.monomials_exponents(n, d; affine::Bool)
#     if affine
#         pred = x -> sum(x) ≤ d
#     else
#         pred = x -> sum(x) == d
#     end
#     E = map(Iterators.filter(pred, Iterators.product(Iterators.repeated(0:d, n)...))) do e
#         collect(e)
#     end

#     sort!(E, lt = _td_order)
#     E
# end
# function td_order(x, y)
#     sx = sum(x)
#     sy = sum(y)
#     sx == sy ? x > y : sx > sy
# end
