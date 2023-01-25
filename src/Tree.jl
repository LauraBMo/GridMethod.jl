

Base.@kwdef mutable struct TreeData{T, d}
    nodes::Vector{Vector{T}} = [zeros(T, d)]
    iter::Int = zero(Int)
end

"""
    Tree{T,d}

Structure representing a tree fractal in `d`-dim space with points of type `T`.
"""
Base.@kwdef struct Tree{T, d}
    data::TreeData{T, d} = TreeData{T, d}()
    ratio::T = T(RATIO[])
end

## Easily create a root Tree.
root(P::Vector{T}, ratio = RATIO[]) where T =
    Tree{T, length(P)}(;
                       data = TreeData{T, length(p)}(; nodes = P),
                       ratio = T(ratio)
                       )
root(::Type{T}, d::Int = 2, ratio = RATIO[]) where T = Tree{T, d}(; ratio = T(ratio))
root(d::Int = 2, ratio = RATIO[]) = root(Float64, d, ratio)

# Access to Tree info
nodes(tree::Tree) = tree.data.nodes
iter(tree::Tree) = tree.data.iter
ratio(tree::Tree) = tree.ratio
factor(tree::Tree) = ratio(tree)^(iter(tree)+one(Int))
dirs(::Tree{T, d}) where {T,d} = pmones(T, Val(d))
branches(tree::Tree) = factor(tree).*dirs(tree)

function Base.copy(tree::Tree{T, d}) where {T,d}
    new_data = TreeData{T, d}(copy(nodes(tree)), copy(iter(tree)))
    return Tree{T, d}(new_data, ratio(tree))
end

function copyempty(tree::Tree{T, d}) where {T,d}
    new_data = TreeData{T, d}(; iter = copy(iter(tree)))
    return Tree{T, d}(new_data, ratio(tree))
end

# function nextnodes(tree::Tree; br = factor(tree).*dirs(tree))
function nextnodes(tree::Tree)
    ## Compute branches once for the whole tree.
    br = branches(tree)
    return [q for p in nodes(tree) for q in [p] .+ br]
end

upiter!(tree::Tree) = tree.data.iter += one(Int)
# upbranches!(tree::Tree) = tree.data.branches = factor(tree).*dirs(tree)

"""
    uptree!(tree, nodes = nextnodes(tree))

Stores `nodes` in `tree.data.nodes` and updates `tree.data.iter`.
Calling `uptree!(tree)` computes next iteration of `tree` and stores it in `tree`.
"""
function uptree!(tree::Tree, nodes = nextnodes(tree))
    @pack! tree.data = nodes
    upiter!(tree)
    return tree
end

"""
    uptree(tree)

Returns the next iteration of `tree`.
"""
function uptree(tree::Tree{T, d}) where {T, d}
    return uptree!(copyempty(tree), nextnodes(tree))
end
nexttree = uptree

function tree(root::Tree, N)
    tree = root
    for _ in 1:N
        uptree!(tree)
    end
    return tree
end

"""
    tree(::Type{T} = Float64, d = 2, N = 0) where T

Returns the `N`th leaves of a `length(dirs(T, d))`-branched tree in a `d`-dimensional space with root the zero point.
"""
tree(::Type{T}, d, N) where T = tree(root(T, d), N)
tree(d = 2, N = 0) = tree(root(d), N)

"""
    findnextleaves(i, tree::Tree)

Returns the `length(dirs(tree))` nodes of `tree` which branched from the `i`-th node of the previous tree in the recursive construction.
Assumes nodes ordered as by [`uptree!`](@ref).
"""
function findnextleaves(i, tree::Tree)
    span = length(dirs(tree))
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
