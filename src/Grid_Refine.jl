# This code is written in Julia and defines several functions that are used to compute the refinement of an initial grid satisfying a condition (e.g., Han condition).

# The grid_point struct represents a point in the grid, with fields x, fx, and step. x is a vector of coordinates, fx is the value of the function at that point, and step is an integer representing the step at which the point was added to the grid. There are also several functions defined for creating and manipulating grid_point objects:

# * grid_point(x, fx, step) creates a new grid_point object with the given coordinates, function value, and step.
# * grid_point(fun::Function, x, step) creates a new grid_point object with the given function, coordinates,
#   and step, where the function value is computed by calling fun(x...).
# * point(g::grid_point) returns the coordinates of the grid_point object g.
# * fpoint(g::grid_point) returns the function value of the grid_point object g.
# * radius(g::grid_point) returns the radius of the grid_point object g, which is computed as 1/2^i,
#   where i is the step at which the point was added to the grid.

# The cube_nthinitgrid function takes a function fun, an integer dim, and an integer n, and returns a list of grid_point objects that form an initial grid with 2^n points in each dimension. The coord_opt argument is a function that is used to compute the coordinates of the points in the grid.

# The cube_pushsubdivided! function takes a list G of grid_point objects, a point p, an integer step, and a function fun, and pushes into G the points obtained by subdividing the point p into 2^dim points.

# The refine_grid function takes a function fun, a parameter C, a parameter m, and an integer dim, and returns a list of grid_point objects that form a refined grid satisfying the Han condition.
# The function works by iterating through the list G and checking each point in turn.
# If the point's function value is less than m, a warning is printed and m is returned.
# If the point satisfies the Han condition (as determined by the isfine function), it is added to the list H. Otherwise, the point is subdivided and the resulting points are added to G.
# This process continues until G is empty.

# The refine_grid! function is similar to refine_grid, but instead of returning a new list of grid_point objects, it modifies the input list G in place.

"""
    grid_point{T}

Structure representing a node in the grid.

# Fields

- `x`: Vector of coordinates of the point.
- `fx`: Value of the function at x.
- `step`: Step of the subdivided point.
"""
struct GridPoint{T, d}
    x::Vector{T}
    fx::T
    step::Int
end

function GridPoint(x, fx, step)
    T = promote_type(eltype(x), typeof(fx))
    return GridPoint{T}(convert(Vector{T}, x), convert(T, fx), Int(step))
end
GridPoint(fun::Function, x, step) = GridPoint(x, fun(x...), step)
# All the evaluations of the function `fun` are done via this function.

step_to_factor(::Type{T}, step) where T = T(^(RATIO[], step)) # 1/2^i
step_to_factor(step) = step_to_factor(Float64, step)

"""
    factor(g::grid_point)

Returns the radius of the `grid_point` object `g`, which is computed as `1/2^i`, where `i` is the step at which the point was added to the grid.
"""
factor(g::GridPoint{T, d}) where {T, d} = step_to_factor(T, g.step)

Base.eltype(::GridPoint{T, d}) where {T, d} = T
function Base.show(io::IO, g::GridPoint)
    print(io, "$(eltype(g))[")
    print(io, join(g.x, ", "))
    print(io, "], $(g.fx), $(g.step)")
end

function grid(tree::Tree{T, d}, fun) where {T, d}
    return [GridPoint(fun, p, iter(tree)+one(Int)) for p in nodes(tree)]
end

"""
    refine_grid(T=Float64, fun, C, m, dim; step₀ = Int(ceil(log2(C))),
                G = cube_nthinitgrid(T, fun, dim, step₀), isfine = _isHan,
                pushsubdivided! = cube_pushsubdivided!)

Returns a list of `grid_point` objects that form a refined grid satisfying the Han condition.

# Arguments:
 - `T` (optional): Numeric type to represent the coordinates of the points.
 - `fun`: The function to be evaluated at each point in the grid.
 - `C`: A parameter used in the Han condition.
 - `m`: A parameter used to determine when to stop refining the grid.
 - `dim`: The dimension of the grid.

# Keyword arguments:
 - `step₀`: The initial step of the grid.
 - `G`: The initial grid.
 - `isfine`: A function that determines whether a point satisfies the Han condition.
 - `pushsubdivide!`: A function that pushes the points obtained by subdividing a region into a list of points.
"""
function refine_grid(::Type{T}, fun, C, m, dim::Int; step₀::Int = Int(ceil(log2(C))),
                     G = grid(tree(T, dim, step₀), fun), isfine = _isHan,
                     pushsubdivided! = cube_pushsubdivided!) where T
    H = eltype(G)[]
    while !(isempty(G)) # That is, while G is not the empty array.
        g = pop!(G) # Removes an element of G and stores it at g.
        @unpack x, fx, step = g
        if fx < m
            @warn "Small norm" fx, m
            return m # WARN Function returning two different types
        end
        println("Step:", step, "($(step_to_factor(step)))", " C:", C, "--", step_to_factor(step) * C, isfine(fx, step, C), fx)
        if isfine(fx, step, C)
            push!(H, g)
        else
            br = branches()
            pushsubdivided!(G, fun, x, br)
        end
    end
    println(step₀)
    return H
end
refine_grid(fun, C, m, dim; step₀ = Int(ceil(log2(C))),
            G = cube_nthinitgrid(fun, dim, step₀), isfine = _isHan,
            pushsubdivided! = cube_pushsubdivided!) =
                refine_grid(Float64, fun, C, m, dim, step₀ = step₀, G = G, isfine = isfine,
                            pushsubdivided! = pushsubdivided!)

function findallandnot(f::Function, A)
    I = findall(f, A)
    return I, setdiff(eachindex(A), I)
end

function push_subdivided!_or_fine!(H, G, br, isfine, m, C, fun)
    g = pop!(H)
    @unpack x, fx, step = g
    if fx < m
        @warn "Small norm" fx, m
        return m # WARN Function returning two different types
    end
    if isfine(fx, step, C)
        push!(G, g)
    else
        pushsubdivided!(H, fun, x, br)
    end
end

checknorms(g, m) = (g.fx < m)
checknorms(G, m) = any(broadcast(checknorms, G, m))

function fine_grid(fun, C, m, dim; isfine = _isHan,
                      pushsubdivide! = cube_pushsubdivided!)
    I = map(isfine, G)
    G, H = G[I], G[(!).(I)]
    step = 1
    while !(isempty(H))
        if checknorms(H, m)
            @warn "Small norm" fx, m
            return m
        end
        br = branches(G, i) # Compute branches for the whole set of pt with step `i`.
        subdivide!(H, G, m, br)
        end
    end
end

function refine_grid!(G, fun, C, m, dim; isfine = _isHan,
                      pushsubdivide! = cube_pushsubdivided!)
    if checknorms(G, m)
        @warn "Small norm" fx, m
        return m
    end
    I = map(isfine, G)
    G, H = G[I], G[(!).(I)]
    while !(isempty(H))
        br = branches(G, H)
        # subdivide!(H, br)
        # if checknorms(H, m)
        #     @warn "Small norm" fx, m
        #     return m
        # end
        # I = map(isfine, H)
        # append!(G, H[I])
        # H = H[J]
        for _ in 1:length(H)
            g = popfirst!(H)
            @unpack x, fx, step = g
            if fx < m
                @warn "Small norm" fx, m
                return m
            end
            if isfine(fx, step, C)
                push!(G, g)
            else
                pushsubdivided!(H, fun, x, br[step])
            end
        end
    end
end

eachstep(H) = unique(map(h -> h.step, H))
branches(G, H) = Dict([(i, branches(G,i)) for i in eachstep(H)])

"""
    refine_grid!(G, fun, C, m, dim; isfine = _isHan,
                 pushsubdivide! = cube_pushsubdivided!)

Modifies the input list `G` in place to form a refined grid satisfying the Han condition.

# Arguments:
 - `G`: A list of `grid_point` objects.
 - `fun`: The function to be evaluated at each point in the grid.
 - `C`: A parameter used in the Han condition.
 - `m`: A parameter used to determine when to stop refining the grid.
 - `dim`: The dimension of the grid.

# Keyword arguments:
 - `isfine`: A function that determines whether a point satisfies the Han condition.
 - `pushsubdivide!`: A function that pushes the points obtained by subdividing a region into a list of points.
"""
function refine_grid!(G, fun, C, m, dim; isfine = _isHan,
                      pushsubdivide! = cube_pushsubdivided!)
    return G = refine_grid(fun, C, m, dim; step₀ = 0, G = G, isfine = isfine,
                           pushsubdivide! = pushsubdivide!)
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
