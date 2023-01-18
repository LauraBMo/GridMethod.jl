
# The Han_min function takes a list H of grid_point objects, a parameter L, and a parameter C, and returns the minimum function value in H (scaled by a factor of γ = 1 - inv(C) * L), as well as the index at which this minimum value was found.

"""
    cube_nthinitgrid(T=Float64, fun, dim, n; ratio = 0.5, root = tree_root(T, dim),
                     dirs = pmones(T, dim))

Returns a list of `grid_point` objects that form an initial grid with `2^n` points in each dimension.

# Arguments:
 - `T` (optional): Numeric type to represent the coordinates of the points.
 - `fun`: The function to be evaluated at each point in the grid.
 - `dim`: The dimension of the grid.
 - `n`: The number of times the grid should be recursively divided.

# Keyword arguments:
 - `ratio`: The ratio between the length of a branch and of the previous.
 - `root`: The root of the grid (a point in `dim`-dimensional space).
 - `coord_opt`: The possible directions in which the next points can be generated.
"""
function cube_nthinitgrid(::Type{T}, fun, dim, n; ratio = 0.5, root = tree_root(T, dim),
                          dirs = pmones(T, dim)) where T
    tree = tree_nthleaves(dim, n; ratio = ratio, root = root, dirs = dirs)
    return [grid_point(fun, p, n) for p in tree]
end
cube_nthinitgrid(fun, dim, n; ratio = 0.5, root = tree_root(dim), dirs = pmones(dim)) =
    cube_nthinitgrid(Float64, fun, dim, n; ratio = ratio, root = root, dirs = dirs)

"""
    cube_pushsubdivided!(G, p, step, fun; dim = length(p),
                         dirs = pmones(eltype(p), dim))

Pushes the points obtained by subdividing the region around `p` into `G`.

# Arguments:
 - `G`: A list of `grid_point` objects, where to push the new points.
 - `p`: A point in `dim`-dimensional space.
 - `step`: The step at which the point was added to the grid.
 - `fun`: The function to be evaluated at each point.

# Keyword arguments:
 - `dim`: The dimension of the grid.
 - `coord_opt`: The possible directions in which the next points can be generated.
"""
function cube_pushsubdivided!(G, p, step, fun; dim = length(p),
                              dirs = pmones(eltype(p), dim))
    ratio = eltype(p)(step_to_ratio(step + 1))
    for q in tree_nextleaves(p, ratio; dim = dim, dirs = dirs)
        push!(G, grid_point(fun, q, step + 1))
    end
end

"""
    _isHan(fx, step, C)

Returns `true` if `C * step_to_ratio(step) < fx`, and `false` otherwise.
"""
_isHan(fx, step, C) = C * step_to_ratio(step) < fx

"""
    Han_min(H, L, C)

Returns the minimum function value in `H` scaled by a factor of `γ`, where `γ = 1 - inv(C) * L`.
"""
function Han_min(H, L, C)
    L < C || @warn "Required condition 0<1-inv(C)L not satisfied." L, C, inv(C) * L
    γ = 1 - inv(C) * L # This is positive because precondition L < C
    minfx, n = findmin(fpoint, H)
    return γ * minfx, n
end

#### Original code below:
#### Explanation and documentation generated with Chat GPT.
#### Steps:
####   1.- "Explain the following code written in Julia language used to compute
####        [explain something, use keywords]
####        ```julia
####        [code goes here]
####        ```"
####   2.- "Use the previous explanation to document de code.
####        Here it is an example of a documented function in Julia language:
####        ```julia
####        [example following your style here]
####        ```
####        "
####   3.- Copy-paste, remove extra sentences, rephrase some sentences... That's pretty much it!
# struct grid_point{T}
#     x::Vector{T}
#     fx::T
#     step::Int
# end
# function grid_point(x, fx, step)
#     U = promote_type(eltype(x), typeof(fx))
#     return grid_point{U}(convert(Vector{U}, x), convert(U, fx), Int(step))
# end
# grid_point(fun::Function, x, step) = grid_point(x, fun(x...), step)

# point(g::grid_point) = g.x
# fpoint(g::grid_point) = g.fx
# radius(g::grid_point) = step_to_ratio(g.step)
# Base.eltype(::grid_point{T}) where {T} = T

# function cube_nthinitgrid(fun, dim, n; ratio = 0.5, root = tree_root(dim),
#                           coord_opt = coord_opt(dim))
#     tree = tree_nthleaves(dim, n; ratio = ratio, root = root, coord_opt = coord_opt)
#     return [grid_point(fun, p, n) for p in tree]
# end

# function cube_pushsubdivided!(G, p, step, fun; dim = length(p),
#                               coord_opt = coord_opt(eltype(p), dim))
#     ratio = step_to_ratio(step + 1)
#     for q in tree_nextleaves(p, ratio; dim = dim, coord_opt = coord_opt)
#         push!(G, grid_point(fun, q, step + 1))
#     end
# end

# _isHan(fx, step, C) = C * step_to_ratio(step) < fx

# function refine_grid(fun, C, m, dim; step₀ = Int(ceil(log2(C))),
#                      G = cube_nthinitgrid(fun, dim, step₀), isfine = _isHan,
#                      pushsubdivided! = cube_pushsubdivided!)
#     println(step₀)
#     H = eltype(G)[]
#     while !(isempty(G)) # That is, while G is not the empty array.
#         g = pop!(G) # Removes an element of G and stores it at g.
#         @unpack x, fx, step = g
#         if fx < m
#             @warn "Small norm" fx, m
#             return m # WARN Function returning two different types
#         end
#         if isfine(fx, step, C)
#             push!(H, g)
#         else
#             pushsubdivided!(G, x, step, fun)
#         end
#     end
#     return H
# end
# function refine_grid!(G, fun, C, m, dim; isfine = _isHan,
#                       pushsubdivide! = cube_pushsubdivided!)
#     return G = refine_grid(fun, C, m, dim; step₀ = 0, G = G, isfine = isfine,
#                            pushsubdivide! = pushsubdivide!)
# end

# # Question: IMPLEMENT finding the min in H's 'while' building loop?
# function Han_min(H, L, C)
#     L < C || @warn "Required condition 0<1-inv(C)L not satisfied." L, C, inv(C) * L
#     γ = 1 - inv(C) * L # This is positive because precondition L < C
#     minfx, n = findmin(fpoint, H)
#     return γ * minfx, n
# end
