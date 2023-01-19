
"""
    C(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
function C(F; kwargs...)
    norm1₀, Δ₀ = norm1(F), Δ(F)
    return C(norm1₀, Δ₀)∘eval_and_J(F; kwargs...)
end

function C(norm1, Δ)
    (xfJ) -> begin
        # unpack
        (x, fx, Jx) = xfJ
        # M = inv(Jx*D(x))*Δ
        M = (Jx*D(x))*inv(Δ)
        # σ = inv(maximum(M))
        σ = maximum(M)
        f2 = maximum(abs.(fx))
        return norm1/max(f2, σ)
    end
end

"""
    K(F; compiled = true, expanded = false)

Returns a function that calculates the condition number of a system of polynomials `F`.

# Keyword arguments:
 - `compiled`: A boolean indicating whether the system should be compiled before being passed to the condition number function.
"""
function K(F; kwargs...)
    Wnorm₀, Δ₀ = Wnorm(F; kwargs...), Δ(F)
    return K(Wnorm₀, Δ₀)∘eval_and_J(F; kwargs...)
end

function K(Wnorm, Δ)
    (xfJ) -> begin
        # unpack
        (x, fx, Jx) = xfJ
        M = inv(Δ)*Jx*D(x)
        σ = last_sval(M, ε)
        f2 = p_norm(fx)^2
        return Wnorm/sqrt(f2 + σ)
    end
end
