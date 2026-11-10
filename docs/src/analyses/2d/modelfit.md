# Custom models

`modelfit2d` fits the amplitudes of each peak across a series to an equation you give.
`x` is the variable that changes from plane to plane, and every other symbol is a
parameter to fit.

## Usage

```julia
using NMRAnalysis
modelfit2d(112:115;
           x=[0.1, 0.2, 0.3, 0.4],
           equation="A*sin(J*x)",
           parameters=["A" => 40.0, "J" => 0.5])
```

`parameters` gives each parameter's starting value. The equation can use `sin`, `cos`,
`exp`, `log`, `sqrt` and arithmetic.

## Output

The fitted parameters, named as in `parameters`, with uncertainties. Plot one against
residue number with `summaryplot("out"; param=:J)`.
