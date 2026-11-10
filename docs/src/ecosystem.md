# Ecosystem

NMRAnalysis.jl is part of a suite of Julia packages for NMR data handling developed
by the [Waudby lab](https://waudbylab.org):

- **[NMRTools.jl](https://waudbylab.org/NMRTools.jl)**: the foundation for
  NMR data import and processing in Julia. NMRAnalysis.jl is built on top of NMRTools
  for all data loading, axis handling, and spectral processing.

- **[NMRExampleData.jl](https://waudbylab.org/NMRExampleData.jl)**: the example datasets
  used throughout this documentation, with saved analyses of many of them. `example(name)`
  downloads a dataset the first time you use it, and `exampleanalysis(name)` gives its
  saved analysis, to open with `reopen`. See the
  [catalogue](https://waudbylab.org/NMRExampleData.jl/stable/catalogue/) for what's
  available.

- **[NMRScreen.jl](https://waudbylab.org/NMRScreen.jl)**: tools for
  fragment and ligand screening by NMR, including automated analysis of large
  compound libraries.
