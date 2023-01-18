
norm1_poly(poly) = sum(abs.(last.(_exposcoeffs(poly))))
norm1(F) = maximum(norm1_poly.(_get_polys(F)))

"""
    C(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
function C(F; kwargs...)
    norm1₀, Δ₀ = norm1(F), Δ(F)
    return C(norm1₀, Δ₀)∘_eval_and_J(F; kwargs...)
end

function C(norm1, Δ)
    (xfJ) -> begin
        # unpack
        (x, fx, Jx) = xfJ
        # M = inv(D(Jx, x))*Δ
        M = (D(Jx, x))*inv(Δ)
        # σ = inv(maximum(M))
        σ = maximum(M)
        f2 = maximum(abs.(fx))
        return norm1/max(f2, σ)
    end
end

"""
    Wnorm_term(I, C)

Returns the term of the Weyl norm corresponding to the exponents `I` and coefficient `C`.
"""
function Wnorm_term(I, C)
    _norm(C)^2*inv(Combinatorics.multinomial(I...))
end
Wnorm_term(IC) = Wnorm_term(IC...)

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
function Wnorm(F; vars = _get_vars(F), expanded = false, kwargs...)
    sum(f -> Wnorm_poly(f; vars = vars, expanded = expanded), _get_polys(F))
end

"""
    K(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
function K(F; kwargs...)
    Wnorm₀, Δ₀ = Wnorm(F; kwargs...), Δ(F)
    return K(Wnorm₀, Δ₀)∘_eval_and_J(F; kwargs...)
end

function K(Wnorm, Δ)
    (xfJ) -> begin
        # unpack
        (x, fx, Jx) = xfJ
        M = inv(Δ)*D(Jx, x)
        σ = last_sval(M, ε)
        f2 = _norm(fx)^2
        return Wnorm/sqrt(f2 + σ)
    end
end
