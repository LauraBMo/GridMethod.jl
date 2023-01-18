## The normsMTK.jl file contains functions for calculating the condition number of a system of polynomials using the ModelingToolkit.jl (MTK) package.

_degree(M::Symbolics.Symbolic, v::AbstractVector) = [Symbolics.degree(M, x) for x in v]
function _exposcoeffs(poly::Symbolics.Num;
                      vars = Symbolics.get_variables(poly),
                      expanded = false,
                      )
    p = expanded ? Symbolics.value(poly) : Symbolics.value(Symbolics.expand(poly))
    return [(_degree(e, vars), c) for (e, c) in p.dict]
end

_get_vars(sys::MTK.NonlinearSystem) = MTK.states(sys)
_get_polys(sys::MTK.NonlinearSystem) = [Symbolics.Num(eq.rhs) for eq in MTK.equations(sys)]
_degrees(sys::MTK.NonlinearSystem) = Symbolics.degree.(_get_polys(sys))

function K(sys::MTK.NonlinearSystem)
    f = MTK.NonlinearFunction(sys;# dvs = states(sys), ps = parameters(sys), u0 = nothing;
                              version = nothing,
                              jac = true,
                              eval_expression = true,
                              sparse = false,
                              simplify = false,
                              )
    eval_and_J(x) = (x, f.f(x,[]), f.jac(x,[]))
    Wnorm₀, Δm1₀ = Wnorm(sys), Δm1(sys)
    return K(Wnorm₀, Δm1₀)∘eval_and_J
end
