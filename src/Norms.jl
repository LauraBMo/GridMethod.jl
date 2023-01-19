
norm1_term(IC) = abs(last(IC))

"""
    norm1_poly(poly; vars, expanded = false)

Returns the sum of the coefficients of a polynomial `poly`.

# Keyword arguments:
 - `vars`: The variables of the polynomial.
 - `expanded`: A boolean indicating whether the polynomial should be expanded before calculating the Weyl norm.
"""
function norm1_poly(poly;
                    vars,
                    expanded = false,
                    )
    sum(norm1_term, exposcoeffs(poly; vars = vars, expanded = expanded))
end

"""
    norm1(F; vars = get_vars(F), expanded = true)

Returns the 1-norm of a system of polynomials `F`.

# Keyword arguments:
 - `vars`: The variables of the system.
 - `expanded`: A boolean indicating whether the polynomials should be expanded before calculating the Weyl norm.
"""
function norm1(F; vars = get_vars(F), expanded = false, kwargs...)
   return maximum(f -> norm1_poly(f; vars = vars, expanded = expanded), get_polys(F))
end

"""
    Wnorm_term(I, C)

Returns the term of the Weyl norm corresponding to the exponents `I` and coefficient `C`.
"""
function Wnorm_term(IC)
    (I, C) = IC
    return p_norm(C)^2*inv(Combinatorics.multinomial(I...))
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
    sum(Wnorm_term, exposcoeffs(poly; vars = vars, expanded = expanded))
end

"""
    Wnorm(F; vars = _get_vars(F), expanded = true)

Returns the Weyl norm of a system of polynomials `F`.

# Keyword arguments:
 - `vars`: The variables of the system.
 - `expanded`: A boolean indicating whether the polynomials should be expanded before calculating the Weyl norm.
"""
function Wnorm(F; vars = get_vars(F), expanded = false, kwargs...)
    sum(f -> Wnorm_poly(f; vars = vars, expanded = expanded), get_polys(F))
end
