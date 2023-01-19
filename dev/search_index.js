var documenterSearchIndex = {"docs":
[{"location":"utils/#Utils","page":"Utils","title":"Utils","text":"","category":"section"},{"location":"norms_const/#Norms-and-constants","page":"Norms and constants","title":"Norms and constants","text":"","category":"section"},{"location":"tree/#Cross-Tree-fractal","page":"Tree fractal","title":"Cross-Tree fractal","text":"","category":"section"},{"location":"tree/","page":"Tree fractal","title":"Tree fractal","text":"The nodes of the grid form the leaves of a cross-formed tree fractal.","category":"page"},{"location":"tree/","page":"Tree fractal","title":"Tree fractal","text":"Modules = [GridMethod]\nPages   = [\"Tree.jl\"]","category":"page"},{"location":"tree/#GridMethod.pmones-Union{Tuple{T}, Tuple{Type{T}, Any}} where T","page":"Tree fractal","title":"GridMethod.pmones","text":"pmones(::Type{T}=Float64, dim)\n\nReturns an iterator which generates all possible dim-dimensional vectors that have all components either one or -one.\n\n\n\n\n\n","category":"method"},{"location":"tree/#GridMethod.tree_nextleaves-Tuple{Any, Any}","page":"Tree fractal","title":"GridMethod.tree_nextleaves","text":"tree_nextleaves(p, ratio; dim = length(p), dirs = pmones(eltype(p), dim))\n\nReturns an iterator which generates the next leaves of the tree starting from p with a length ratio of r.\n\nKeyword arguments:\n\ndim: The dimension of the tree.\ndirs: The directions in which the next leaves will be generated.\n\n\n\n\n\n","category":"method"},{"location":"tree/#GridMethod.tree_nthleaves-Tuple{Any, Any}","page":"Tree fractal","title":"GridMethod.tree_nthleaves","text":"tree_nthleaves(dim, n; ratio = 0.5, root = tree_root(dim),\n               dirs = pmones(eltype(p), dim))\n\nComputes the nth leaves of a 2^dim-branched tree in a dim-dimensional space, where dim is the dimension of the tree (i.e. 2 for a 2D tree, 3 for a 3D tree), and n is the number of times the tree should be recursively branched.\n\nKeyword arguments:\n\nratio: The ratio between the length of a branch and of the previous.\nroot: The root of the tree (a point in dim-dimensional space).\ndirs: The possible directions in which the next leaves can be generated.\n\n\n\n\n\n","category":"method"},{"location":"tree/#GridMethod.tree_root-Union{Tuple{T}, Tuple{Type{T}, Any}} where T","page":"Tree fractal","title":"GridMethod.tree_root","text":"tree_root(::Type{T}=Float64, dim)\n\nReturns an array containing a single zero-vector of length dim and type T. It is meant to be the root node of a tree fractal.\n\n\n\n\n\n","category":"method"},{"location":"grid_refine/#The-grid-and-refining-it","page":"Grid and refine","title":"The grid and refining it","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = GridMethod","category":"page"},{"location":"#GridMethod","page":"Introduction","title":"GridMethod","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Documentation for GridMethod.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This is a project aiming at implementing the Grid Method theoretically developed by Bürgisser, Cucker, Krick, Lairez, Shub and Tonelli-Cueto.  The aim of this project is to develop parallelizable algorithms for solving real polynomial systems and computing the topology of real algebraic and semialgebraic sets...","category":"page"},{"location":"#Contents","page":"Introduction","title":"Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Problem formulation","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\n    \"modelingkits.md\",\n]\nDepth = 2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The grid","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\n    \"grid_refine.md\",\n    \"cubeHan.md\",\n]\nDepth = 2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Cross-tree fractal","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\n    \"tree.md\",\n]\nDepth = 2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Norms and constants","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\n    \"norms_const.md\",\n]\nDepth = 2","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Utils","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\n    \"utils.md\",\n]\nDepth = 2","category":"page"},{"location":"modelingkits/#Modeling-kits","page":"Modeling Kits","title":"Modeling kits","text":"","category":"section"},{"location":"modelingkits/","page":"Modeling Kits","title":"Modeling Kits","text":"We do not implement a modeling language. Instead, we offer easily to build compatibility with your favorite system. Available packages are HomotopyContinuation.ModelKit.jl and ModelingToolkit.jl. See their config files HC.ModelKit.jl and ModelingToolkit.jl.","category":"page"},{"location":"modelingkits/","page":"Modeling Kits","title":"Modeling Kits","text":"Modules = [GridMethod]\nPages   = [\"HC.ModelKit.jl\", \"ModelingToolkit.jl\"]","category":"page"},{"location":"cubeHan/#Dividing-the-cube-under-Han-condition","page":"Cube – Han","title":"Dividing the cube under Han condition","text":"","category":"section"}]
}
