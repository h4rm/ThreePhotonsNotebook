IJulia Notebook accompanying the following paper:

You can install IJulia, PyPlot and seaborn with:

```julia
ENV["JUPYTER"] = ""
ENV["PYTHON"] = ""
Pkg.add("IJulia")
Pkg.add("PyPlot")
import Conda
Conda.add("seaborn")
```

You will further have to install the [ThreePhotons.jl](https://github.com/h4rm/ThreePhotons.jl) package.
