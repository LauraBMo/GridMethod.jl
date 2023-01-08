
# The tree_root function creates a list containing a single element, which is an array of zeros with dim dimensions. This is meant to be used as the root node of a tree fractal.
# The coord_opt function returns a tuple of pairs, with dim pairs in total. Each pair consists of one and -one of the specified type. This is meant to be used to define the possible directions in which the branches of the tree fractal can extend.
# The tree_nextleaves function takes in a point p and a ratio ratio, and returns an iterator that generates all possible next points for the tree fractal, based on the coord_opt tuple.
# The tree_nthleaves function generates the nth level leaves of a tree fractal with dim dimensions. It generates new points for the fractal iteratively.

"""
    tree_root(::Type{T}, dim) where {T} = [zeros(T, dim)]

Returns an array containing a single zero-vector of length `dim` and type `T`, which defaults to `Float64`.
It is meant to be used as the root node of a tree fractal.
"""
tree_root(::Type{T}, dim) where {T} = [zeros(T, dim)]
tree_root(dim) = tree_root(Float64, dim)

"""
    coord_opt(::Type{T}, dim) where {T} = ntuple(_ -> (one(T), -one(T)), dim)

Returns an iterator which generates all possible `dim`-dimensional vectors that have all components either one or -one.
The `eltype` of the vector is `T`, which defaults to `Float64`.
"""
coord_opt(::Type{T}, dim) where {T} = ntuple(_ -> (one(T), -one(T)), dim)
coord_opt(dim) = coord_opt(Float64, dim)

"""
    tree_nextleaves(p, ratio; dim = length(p), coord_opt = coord_opt(eltype(p), dim))

Returns an iterator which generates the next leaves of the tree starting from `p` with a length ratio of `r`.

# Keyword arguments:
 - `dim`: The dimension of the tree.
 - `coord_opt`: The possible directions in which the next leaves can be generated.
"""
function tree_nextleaves(p, ratio; dim = length(p), coord_opt = coord_opt(eltype(p), dim))
    # Generate an iterator of possible directions using coord_opt
    dir_vec = Iterators.ProductIterator(coord_opt)
    # Define a function that takes a direction vector and returns the next leave in that direction
    f = v -> p + ratio * collect(v)
    # Return an iterator that applies the function f to the direction vectors
    return Iterators.map(f, dir_vec)
end

"""
    tree_nthleaves(dim, n; ratio = 0.5, root = tree_root(dim),
                   coord_opt = coord_opt(dim))

Computes the `n`th leaves of a 2^`dim`-branched tree in a `dim`-dimensional space, where `dim` is the dimension of the tree (i.e. 2 for a 2D tree, 3 for a 3D tree), and `n` is the number of times the tree should be recursively branched.

# Keyword arguments:
 - `ratio`: The ratio between the length of a branch and of the previous.
 - `root`: The root of the tree (a point in `dim`-dimensional space).
 - `coord_opt`: The possible directions in which the next leaves can be generated.
"""
function tree_nthleaves(dim, n; ratio = 0.5, root = tree_root(dim),
                        coord_opt = coord_opt(dim))
    tree = root
    r = ratio
    for _ in 1:n
        new_tree = eltype(root)[]
        for p in tree
            # Get an iterator of the next leaves of the tree starting from p with a length ratio of r
            for q in tree_nextleaves(p, r; dim = dim, coord_opt = coord_opt)
                push!(new_tree, q)
            end
        end
        tree = new_tree
        r *= ratio
    end
    return tree
end
