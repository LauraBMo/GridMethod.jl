module GridMethod

#Indicate modules for use
using UnPack: @unpack
# using LinearAlgebra: dot
# using Parameters: @unpack
# using DocStringExtensions: SIGNATURES, TYPEDEF
# using Requires # for macro @require
# using RecipesBase # For Plots recipes

using Combinatorics
import LinearAlgebra as LA

#Include files with extracode
export rand_poly, Wnorm_poly, Wnorm, K
include("utils.jl")

include("HanSubdivision.jl")
include("cube_tree.jl")


# import HomotopyContinuation as HC
import HomotopyContinuation.ModelKit as MK
include("normsHC.jl")


# include("normsSym.jl")

import Symbolics
import ModelingToolkit as MTK
include("normsMTK.jl")



end
