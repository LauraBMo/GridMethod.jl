
"""
    Tree{T,d}

Structure representing a tree fractal in `d`-dim space with points of type `T`.
"""
Base.@kwdef mutable struct Tree{T,d}
    nodes::Vector{Vector{T}} = [zeros(T, d)]
    iter::Int = 0
end

Tree(::Type{T} = Float64, d = 2) where T = Tree{T, d}()
nodes(tree::Tree) = tree.nodes
Base.length(tree::Tree) = length(tree.nodes)
Base.eltype(::Tree{T, d}) where {T, d} = Tuple{Int, Vector{T}}

function Base.iterate(tree::Tree, state = 1)
    if state > length(tree)
        return nothing
    end
    i = state
    return (i, nodes(tree)[i]), state + 1
end

"""
    pmones(d::Int = 2, ::Type{T} = Float64)

Returns an iterator which generates all possible `d`-dimensional vectors that have all components either one or -one.
"""
pmones(d::Int = 2, ::Type{T} = Float64) where {T} = Iterators.ProductIterator(ntuple(_ -> (one(T), -one(T)), d))

"""
    tree(root::Tree{T, d}, N; ratio = 0.5, dirs = pmones(T, dim))

Returns the `N`th leaves of a `2^d`-branched tree in a `d`-dimensional space with root `root`, where `N` is the number of times the root-tree should be recursively branched into directions `dirs`.

# Keyword arguments:
 - `ratio`: The ratio between the length of the recursively branched.
 - `dirs`: The possible directions in which the next leaves can be generated.
"""
function tree(root::Tree{T, d}, N; ratio = .5, dirs = pmones(T, d)) where {T,d}
    tree, r = root, ratio
    for _ in 1:N
        tree = tree_nexttree(tree, r, dirs)
        r *= ratio
    end
    return tree
end
tree(d::Integer, N; ratio = .5, dirs = pmones(d)) =
    tree(Tree(d), N; ratio = ratio, dirs = dirs)

"""
    tree_nexttree(tree, ratio, dim, dirs = pmones(T, d))

Returns the next iteration of `tree` with a ratio of length `ratio`.
Where `dim` is the dimension of the tree, and `dirs` the directions generating next leaves.
"""
function tree_nexttree(tree::Tree{T, d}, ratio, dirs = pmones(T, d)) where {T,d}
    return tree_nexttree(tree, tree_nextleaves(ratio, dirs))
end
tree_nexttree(tree::Tree{T, d}, nextleaves::Function) where {T, d} =
    Tree{T, d}([q for p in nodes(tree) for q in nextleaves(p)],
         tree.iter + 1
         )

"""
    tree_nextleaves(ratio, dirs)

Returns a function whose input is a point `p` and the output is an iterator which generates the next leaves of the tree starting from `p` with a ratio of length `ratio` into the directions in the iterator `dirs`.
"""
function tree_nextleaves(ratio, dirs)
    # Define a function that takes a direction vector and returns the next leave in that direction
    # f = v -> p + ratio * collect(v)
    # Return an iterator that applies the function f to the direction vectors
    return p -> Iterators.map(v -> p + ratio*collect(v), dirs)
end

"""
    tree_findnextleaves(i, tree::Tree{T, d}) where {T, d}

Returns the `2^d-1` nodes of `tree` which branched from the `i`-th node of the previous tree in the recursive construction.
"""
function tree_findnextleaves(i, tree::Tree{T, d}) where {T, d}
    span = 2^d
    init = span*(i-1)
    return nodes(tree)[init+1:init+span]
end

## For a more granular interaction while building the tree:
# """
#     tree_firstnextleaves!(tree::Tree{T, d}, nextleaves) where {T,d}

# Deletes first node of the tree and inserts its division at the end.
# """
# function tree_firstnextleaves!(tree::Tree{T, d}, nextleaves) where {T,d}
#     p = popfirst!(tree.nodes)
#     for q in nextleaves(p)
#         push!(tree.nodes, q)
#     end
# end

# function tree_newleaves(tree::Tree{T, d}) where {T, d}
#     return tree.nodes[end-(2^d):end]
# end
