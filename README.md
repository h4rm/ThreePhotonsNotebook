Supplementary IJulia Notebook for the paper "*Structure determination from single molecule X-ray scattering with three photons per image*".

In addition to the dependencies below, you will have to install the [https://github.com/h4rm/ThreePhotons.jl](https://github.com/h4rm/ThreePhotons.jl) package (not listed in the Official Julia Repository).

The data that is analyzed in this notebook is not yet released.

Dependencies
=============

[Julia (5.1)](https://julialang.org/downloads/oldreleases.html)

Packages
---------

* [IJulia.jl](https://github.com/JuliaLang/IJulia.jl)
* [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl)
* [LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl)
* [Seaborn](https://seaborn.pydata.org/) (Matplotlib Package provided via Conda)

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
=====================

You can run the notebook by

```julia
using IJulia
notebook() #A browser window should pop up
```
