```@setup treefractal
using Plots; gr()
Plots.reset_defaults()
```

# Cross-Tree fractal in cube

  The nodes of the grid form the leaves of a cross-form tree fractal.
  The `tree` function generates the Nth level leaves of a tree fractal with dim dimensions. 
  It generates new points for the fractal iteratively.
  
#### Tree fractal

  ![Tree fractal in 2D](./assets/tree_lines.gif)

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

function plot_lines!(P, N, r)
    HH = [GridMethod.tree_nthleaves(2, i) for i in 0:N]
    r₀ = r
    @gif for H in HH
        for p in H
            for q in GridMethod.tree_nextleaves(p, r₀)
                pp = [p, q, [NaN, NaN]]
                plot!(P, first.(pp), last.(pp), line = (:black, 1, 0.6))
            end
        end
        r₀ *= r
    end every 1
end

P = plot(1,
         legend = false,
         color = :black,
         xlim = (-1, 1),
         ylim = (-1, 1),
         );
plot_lines!(P, 4, .5)
```
