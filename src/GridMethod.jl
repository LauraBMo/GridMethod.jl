module GridMethod

#Indicate modules for use
using UnPack # @unpack, @pack!
# using DocStringExtensions: SIGNATURES, TYPEDEF
# using RecipesBase # For Plots recipes
# using InvertedIndices

using Combinatorics
import LinearAlgebra as LA

#Include files with extracode
export rand_poly
include("Utils.jl")

export norm1, C, Wnorm, K
include("Norms.jl")
include("Const.jl")

const RATIO = Ref{Float64}(.5)

"""
    set_default_ratio(ratio)

Set default ratio for Tree structure.
"""
function set_default_ratio(ratio)
    RATIO[] = ratio
end

export tree, root, uptree!, uptree, findnextleaves
include("Tree.jl")

include("Grid_Refine.jl")
include("CubeHan.jl")



using Requires

function __init__()
    @require HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327" begin
        # import .HomotopyContinuation as HC
        import .HomotopyContinuation.ModelKit as MK
        include("HC.ModelKit.jl")
    end
    @require ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78" begin
        @eval import Symbolics
        import .ModelingToolkit as MTK
        include("ModelingToolkit.jl")
    end
end


end
