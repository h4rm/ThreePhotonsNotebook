Supplementary IJulia Notebook for the paper "*Structure determination from single molecule X-ray scattering with three photons per image*".

In addition to the dependencies below, you will have to install the [https://github.com/h4rm/ThreePhotons.jl](https://github.com/h4rm/ThreePhotons.jl) package (not listed in the Official Julia Repository).

The Crambin data that is analyzed in this notebook is not yet released, but can be reproduced with the job scripts under `jobs/run_generate_histograms.jl` and `jobs/run_parallel_determination.jl`.

Content
=============

* *Paper.ipynb* - Main IJulia notebook with many of the necessary data processing steps. A lot of the code for postprocessing has been put into `paper.jl` (in the directory) and `data_processing.jl` in the ThreePhoton.jl package.
* *Paper_revision.ipynb* - Data processing discussed in the revision of the paper, including experimental data (Coliphage_PR772), multi-particle processing and Ewald sphere approximations.
* *Beamstop.ipynm* - Calculations concering the beamstop, both in synthetic experiments and real experimental images.

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

You can run the notebook by executing

```julia
using IJulia
notebook() #A browser window should pop up
```
