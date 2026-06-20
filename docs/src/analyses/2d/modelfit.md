# Custom Model Fitting

```@docs
modelfit2d
```

The `modelfit2d` function fits peak amplitudes across a series of 2D spectra to an
arbitrary user-supplied equation. This is useful when none of the built-in models
(exponential decay, recovery) describe your experiment.

The equation is specified as a string using standard mathematical syntax. The variable `x`
represents the independent variable (e.g. a delay time, concentration, or frequency). All
other symbols in the equation are treated as free parameters to be fitted.

<!-- screenshot: assets/modelfit2d.png -->

## Usage

```julia
using NMRAnalysis

# J-modulation experiment: amplitudes oscillate as A·sin(J·x)
modelfit2d(
    ["112/pdata/1", "113/pdata/1", "114/pdata/1", "115/pdata/1"],
    [0.1, 0.2, 0.3, 0.4],        # x values (e.g. delay times in seconds)
    "A*sin(J*x)",                 # model equation
    ["A" => 40.0, "J" => 0.5]    # parameter names and initial values
)
```

Parameters are specified as a vector of `"name" => initial_value` pairs. Good initial
values improve convergence, especially for nonlinear models.

Common mathematical functions available in equations include `sin`, `cos`, `exp`, `log`,
`sqrt`, and standard arithmetic operators.

## Output

Clicking **Save to folder** writes fitted parameter values to `fit-results.txt`. The
column names in the output file match the parameter names given in the `parameters`
argument.
