## The utils.jl file contains utility functions for working with polynomials.

"""
    rand_poly(::Type{T}=Float64, X, d; coeffs = I -> randn(T), n = length(X))

Returns a random polynomial with coefficients of type `T` in variables `X` with degree `d`.

# Keyword arguments:
 - `coeffs`: A function that takes a tuple of exponents and returns a coefficient.
 - `n`: The number of variables.
"""
function rand_poly(::Type{T}, X, d; coeffs = I -> randn(T), n = length(X)) where T
    sum(I -> coeffs(I)*_monomial(X, I), Combinatorics.multiexponents(n, d))
end
rand_poly(X, d; coeffs = I -> randn(Float64), n = length(X)) =
    rand_poly(Float64, X, d; coeffs = coeffs, n = n)

const ε = 1e-6 # bound for meaning zero
function last_sval(M, ε)
    U, svals, Vt = LA.svd(M; alg = QRIteration())
    # svals = LA.svdvals(M; alg = QRIteration())
    return svals[findlast(>(ε), svals)]
end

## From HC.jl
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
    Id(::Type{T}=Float64, n)

Returns the identity matrix of size `n` and type `T`.
"""
function Id(::Type{T}, n) where T
    LA.Diagonal(ones(T, n))
end
Id(n) = Id(Float64, n)
# diagonal(v) = cat(v...; dims = (1, 2)) # Without LinearAlgebra.jl

"""
    _norm(v, p = 2)

Returns the p-norm of a vector `v`.
"""
function _norm(v, p = 2)
    sum(x -> x^p, v)^(1/p)
end

"""
    Δ(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function Δ(F)
    LA.Diagonal(_degrees(F))
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

## For docs
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
