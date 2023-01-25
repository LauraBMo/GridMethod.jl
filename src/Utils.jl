## The utils.jl file contains utility functions for working with polynomials.

"""
    rand_poly(::Type{T}=Float64, X, d; coeffs = I -> randn(T), n = length(X))

Returns a random polynomial of coefficient-type `T` in variables `X` with degree `d`.

# Keyword arguments:
 - `coeffs`: A function that takes a tuple of exponents and returns a coefficient.
 - `n`: The number of variables.
"""
function rand_poly(::Type{T}, X, d; coeffs = I -> randn(T), n = length(X)) where T
    sum(I -> coeffs(I)*monomial(X, I), Combinatorics.multiexponents(n, d))
end
rand_poly(X, d; coeffs = I -> randn(Float64), n = length(X)) =
    rand_poly(Float64, X, d; coeffs = coeffs, n = n)

const ε = 1e-6 # bound for meaning zero
function last_sval(M, ε)
    U, svals, Vt = LA.svd(M; alg = LA.QRIteration())
    # svals = LA.svdvals(M)
    return svals[findlast(>(ε), svals)]
end

## From HC.jl
"""
    monomial(X, I, n = length(X))

Returns the monomial of variables `X` with exponents `I`.

# Keyword arguments:
 - `n`: The number of variables.
"""
function monomial(X, I, n=length(X))
    prod(i -> X[i]^I[i], 1:n)
end

"""
    Id(::Type{T}=Float64, n)

Returns the identity matrix of size `n` and type `T`.
"""
function Id(::Type{T}, n) where T
    LA.Diagonal(ones(T, n))
end
Id(v::AbstractVector) = Id(eltype(v), length(v))
# diagonal(v) = cat(v...; dims = (1, 2)) # Without LinearAlgebra.jl

"""
    p_norm(v, p = 2)

Returns the p-norm of a vector `v`.
"""
function p_norm(v, p = 2)
    sum(x -> x^p, v)^(1/p)
end

"""
    Δ(F)

Returns the degree matrix of a system of polynomials `F`.
"""
function Δ(F)
    LA.Diagonal(degrees(F))
end

"""
    D(x)

Returns the matrix difference between the identity matrix and the outer product of `x` and its conjugate transpose.

# Keyword arguments:
 - `n`: The size of the identity matrix.
"""
function D(x)
    Id(x) - x*transpose(conj(x))
end

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
