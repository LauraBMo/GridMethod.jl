module GridMethod

#Indicate modules for use
using UnPack: @unpack
# using DocStringExtensions: SIGNATURES, TYPEDEF
# using RecipesBase # For Plots recipes

using Combinatorics
import LinearAlgebra as LA

#Include files with extracode
export rand_poly
include("Utils.jl")

export norm1, C, Wnorm, K
include("Norms.jl")
include("Const.jl")

export tree_root, pmones, tree

export Tree, pmones, tree, tree_nexttree, tree_findnextleaves
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
