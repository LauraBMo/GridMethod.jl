
mutable struct TreeData{T, d}
    nodes::Vector{Vector{T}}
    iter::Int
end

"""
    Tree{T,d}

Structure representing a tree fractal in `d`-dim space with points of type `T`.
"""
Base.@kwdef struct Tree{T, d}
    data::TreeData{T, d} = TreeData{T, d}([zeros(T, d)], zero(Int))
    ratio::T = inv(one(T)+one(T))
end

## Easily create a root Tree.
Tree(::Type{T}, d) where T = Tree{T, d}()
Tree(d::Int = 2) = Tree(Float64, d)

# Access to Tree info
nodes(tree::Tree) = tree.data.nodes
iter(tree::Tree) = tree.data.iter
ratio(tree::Tree) = tree.ratio

function Base.copy(tree::Tree{T, d}) where {T,d}
    new_data = TreeData{T, d}(copy(nodes(tree)), copy(iter(tree)))
    return Tree{T, d}(; data = new_data, ratio = ratio(tree))
end

factor(tree::Tree) = ratio(tree)^(iter(tree)+one(Int))

"""
    nextleaves(tree::Tree)

Returns an iterator which generates the directions for the next leaves of the tree.
"""
function nextleaves(tree::Tree{T, d}) where {T, d}
    dirs = Iterators.ProductIterator(ntuple(_ -> (one(T), -one(T)), d))
    return Iterators.map(v -> factor(tree)*collect(v), dirs)
end

# Iterate over a Tree `tree`: `(i, p)` where `p` is the `i`-th node of `tree`.
Base.length(tree::Tree) = length(nodes(tree))
Base.eltype(::Tree{T, d}) where {T, d} = Tuple{Int, Vector{T}}
function Base.iterate(tree::Tree, state = one(Int))
    if state > length(tree)
        return nothing
    end
    i = state
    return (i, nodes(tree)[i]), state + one(Int)
end

# Update (compute next iteration) data in a tree:
upiter!(tree::Tree) = tree.data.iter += one(Int)
setnodes!(tree::Tree, nodes) = tree.data.nodes = nodes

"""
    uptree!(tree)

Returns the next iteration of `tree`.
"""
function uptree!(tree::Tree)
    setnodes!(tree, [q for p in nodes(tree) for q in [p] .+ nextleaves(tree)])
    upiter!(tree)
    return tree
end

"""
    tree(root::Tree, N)

Returns the `N`th leaves of a `2^d`-branched tree in a `d`-dimensional space with root `root`, where `N` is the number of times the root-tree should be recursively branched into directions `dirs`.

# Keyword arguments:
 - `ratio`: The ratio between the length of the recursively branched.
 - `dirs`: The possible directions in which the next leaves can be generated.
"""
function tree(root::Tree, N)
    tree = root
    for _ in 1:N
        uptree!(tree)
    end
    return tree
end
tree(d::Int = 2, N = 0) = tree(Tree(d), N)

"""
    findnextleaves(i, tree::Tree)

Returns the `2^d` nodes of `tree` which branched from the `i`-th node of the previous tree in the recursive construction. Assumes nodes ordered as given by [`uptree!`](@ref).
"""
function findnextleaves(i, tree::Tree)
    span = length(nextleaves(tree))
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
