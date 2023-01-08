## The normsHC.jl file contains functions for calculating the condition number of a system of polynomials using the ModelKit.jl (MK) package.

function _exposcoeffs(poly::MK.Expression;
                      vars = MK.variables(poly),
                      expanded = false,
                      )
    p = expanded ? poly : MK.expand(poly)
    return [(e, MK.to_number(c)) for (e, c) in MK.to_dict(p, vars)]
end

_get_vars(F::MK.System) = MK.variables(F)
_get_polys(F::MK.System) = MK.expressions(F)
_degrees(F::MK.System) = MK.degrees(F)

# K(system, vars) = x -> K(x, K_prepare(system, vars))
function K(F::MK.System; compiled = true)
    f = compiled ? MK.CompiledSystem(F; optimizations = true) : MK.InterpretedSystem(F; optimizations = true)
    Wnorm₀, Δm1₀ = Wnorm(F), Δm1(F)
    function eval_and_J(x)
        u = Vector{Any}(undef, size(F, 1))
        U = Matrix{Any}(undef, size(F))
        MK.evaluate_and_jacobian!(u, U, f, x, nothing)
        Jx = MK.to_smallest_eltype(U)
        fx = MK.to_smallest_eltype(u)
        return (x, fx, Jx)
    end
    return K(Wnorm₀, Δm1₀)∘eval_and_J
end
