
"""
    tree_root(::Type{T}=Float64, dim)

Returns an array containing a single zero-vector of length `dim` and type `T`.
It is meant to be the root node of a tree fractal.
"""
tree_root(::Type{T}, dim) where {T} = [zeros(T, dim)]
tree_root(dim) = tree_root(Float64, dim)

"""
    pmones(::Type{T}=Float64, dim)

Returns an iterator which generates all possible `dim`-dimensional vectors that have all components either one or -one.
"""
pmones(::Type{T}, dim) where {T} = Iterators.ProductIterator(ntuple(_ -> (one(T), -one(T)), dim))
pmones(dim) = pmones(Float64, dim)

"""
    tree_nextleaves(p, ratio; dim = length(p), dirs = pmones(eltype(p), dim))

Returns an iterator which generates the next leaves of the tree starting from `p` with a length ratio of `ratio`.

# Keyword arguments:
 - `dim`: The dimension of the tree.
 - `dirs`: The directions in which the next leaves will be generated.
"""
function tree_nextleaves(p, ratio; dim = length(p), dirs = pmones(eltype(p), dim))
    # Define a function that takes a direction vector and returns the next leave in that direction
    f = v -> p + ratio * collect(v)
    # Return an iterator that applies the function f to the direction vectors
    return Iterators.map(f, dirs)
end

"""
    tree_nexttree(tree, ratio, dim, dirs)

Returns the next iteration of `tree` with a ratio of length `ratio`.
Where `dim` is the dimension of the tree, and `dirs` the directions generating next leaves.
"""
function tree_nexttree(tree, ratio, dim, dirs)
    return [q for p in tree for q in tree_nextleaves(p, ratio; dim = dim, dirs = dirs)]
end

"""
    tree(::Type{T}=Float64, dim, N; ratio = 0.5, root = tree_root(T, dim),
         dirs = pmones(T, dim))

Computes the `N`th leaves of a 2^`dim`-branched tree in a `dim`-dimensional space, where `dim` is the dimension of the tree (i.e. 2 for a 2D tree, 3 for a 3D tree), and `n` is the number of times the tree should be recursively branched.

# Keyword arguments:
 - `ratio`: The ratio between the length of a branch and of the previous.
 - `root`: The root of the tree (a point in `dim`-dimensional space).
 - `dirs`: The possible directions in which the next leaves can be generated.
"""
function tree(::Type{T}, dim, N; ratio = 0.5, root = tree_root(T, dim), dirs = pmones(T, dim)) where T
    tree = root
    r = ratio
    for _ in 1:N
        tree = tree_nexttree(tree, r, dim, dirs)
        r *= ratio
    end
    return tree
end
tree(dim, N; ratio = 0.5, root = tree_root(dim), dirs = pmones(dim)) =
    tree(Float64, dim, N; ratio = ratio, root = root, dirs = dirs)
