IJulia Notebook accompanying the following paper:

In addition to the dependencies below, you will have to install the [ThreePhotons.jl](https://github.com/h4rm/ThreePhotons.jl) package (not listed in the Official Julia Repository).

Dependencies
############

* IJulia.jl
* PyPlot.jl
* LaTeXStrings.jl
* Seaborn (Matplotlib Package provided via Conda)

You can install the dependencies via:

```julia
  ENV["JUPYTER"] = ""
  ENV["PYTHON"] = ""
  Pkg.add("IJulia")
  Pkg.add("PyPlot")
  import Conda
  Conda.add("seaborn")
  Pkg.add("LaTeXStrings")
```

Running the Notebook
#####################

You can run the notebook by

```julia
using IJulia
notebook() #A browser window should pop up
```
