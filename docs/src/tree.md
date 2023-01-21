```@setup treefractal
using Plots; gr()
Plots.reset_defaults()
```

# Cross-Tree fractal in cube

  The nodes of the grid form the leaves of a cross-form tree fractal.
  The `tree` function generates the Nth level leaves of a tree fractal with dim dimensions. 
  It generates new points for the fractal iteratively.
  
#### Tree fractal

  ![Tree fractal in 2D](./assets/tree_lines3.gif)

  *See [Image code](@ref)*

## Documentation

```@autodocs
Modules = [GridMethod]
Pages   = ["Tree.jl"]
Private = false
```

## Image code

Code used to generate image [Tree fractal](@ref).

```@example treefractal
using Plots
using GridMethod

function plot_lines!(P, N, w, dw, c, dc)
    w₀, c₀ = w, c
    tree = GridMethod.tree()
    @gif for _ in 1:N
        old_tree = copy(tree)
        uptree!(tree)
        for (i, p) in old_tree
            for q in findnextleaves(i, tree)
                pp = [p, q, [NaN, NaN]]
                plot!(P, first.(pp), last.(pp), line = (:black, w₀, c₀))
            end
        end
        w₀ += dw
        c₀ += dc
    end every 1
end

P = plot(1,
         legend = false,
         color = :black,
         xlim = (-1, 1),
         ylim = (-1, 1),
         );
plot_lines!(P, 7, .5, 2.3, -.4, .75, -.11)
display(P)
```
