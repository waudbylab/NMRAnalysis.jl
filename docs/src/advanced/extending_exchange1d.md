# Extending Exchange1D

This guide explains how to add new exchange models and experiment types to the Exchange1D module. The module uses Julia's type dispatch system, making it straightforward to extend without modifying existing code.

## Architecture Overview

The Exchange1D module is organised around two abstract type hierarchies:

```
AbstractModel                    AbstractExperiment
├── NoExchangeModel              ├── R1Experiment
├── TwoStateModel                └── CESTExperiment
├── ThreeStateModel
└── TwoStateBindingModel
```

An `ExchangeProblem` combines a vector of experiments with a model, and the `fit()` function performs joint least-squares optimisation across all experiments.

### Module Structure

```
src/exchange1d/
├── Exchange1D.jl          # Module entry point, analysis registration
├── types.jl               # AbstractModel, AbstractExperiment, ExchangeProblem, FitResult
├── models.jl              # Includes all model-*.jl files
├── model-noexchange.jl    # NoExchangeModel implementation
├── model-twostate.jl      # TwoStateModel implementation
├── model-threestate.jl    # ThreeStateModel implementation
├── model-twostatebinding.jl # TwoStateBindingModel implementation
├── experiments.jl         # load_experiment(), field_label(), includes expt-*.jl
├── expt-r1.jl             # R1Experiment implementation
├── expt-cest.jl           # CESTExperiment implementation
├── params.jl              # Parameter assembly (defaultparams for ExchangeProblem)
├── problem.jl             # simulate!, residuals, fit(), plot_result
├── liouvillian.jl         # Bloch-McConnell Liouvillian construction
├── interface.jl           # Interactive text-based UI
├── results.jl             # FitResult display
├── fitting-with-errors.jl # Uncertainty propagation
└── misc.jl                # Output folder utilities
```

### Parameter System

Parameters are stored in a `ComponentArray` with three sections:

```julia
params = ComponentArray(;
    model    = ...,   # Exchange parameters (kex, pB, Kd, etc.)
    spin     = ...,   # Chemical shifts and relaxation rates
    nuisance = ...,   # Per-experiment amplitude/correction factors
)
```

Field-dependent parameters are tagged with the magnetic field strength using `field_label()`, which converts e.g. `14.1` to `:14p1T`. This ensures experiments at different fields get separate relaxation rate parameters while experiments at the same field share them.

## Adding a New Model

New models are automatically discovered at runtime via Julia's type system — no registration step is needed. Any concrete subtype of `AbstractModel` will appear in the model selection menu.

### Step 1: Create the model file

Create a new file `src/exchange1d/model-yourmodel.jl`.

### Step 2: Define the struct

```julia
struct YourModel <: AbstractModel
end
```

For models involving multiple molecules (e.g. binding), add a `moleculemap` field:

```julia
struct YourModel <: AbstractModel
    moleculemap::Dict{Symbol,String}
end
YourModel() = YourModel(Dict{Symbol,String}())
```

### Step 3: Implement required functions

Every model must implement the following eight functions:

#### `modelname(::YourModel) -> String`

Display name shown in the model selection menu.

```julia
modelname(::YourModel) = "Your exchange model"
```

#### `nstates(::YourModel) -> Int`

Number of exchanging states.

```julia
nstates(::YourModel) = 2
```

#### `states(::YourModel) -> Vector{String}`

Labels for each state, used in parameter display. Length must equal `nstates`.

```julia
states(::YourModel) = ["A", "B"]
```

#### `nmolecules(::YourModel) -> Int`

Number of distinct molecular species. Use `1` for intramolecular exchange, `2` for bimolecular binding.

```julia
nmolecules(::YourModel) = 1
```

#### `molecules(::YourModel) -> Dict{Symbol,String}`

Maps role symbols to human-readable descriptions. Used to prompt the user for molecule assignment when `nmolecules > 1`.

```julia
molecules(::YourModel) = Dict(:A => "observed")
```

#### `defaultparams(::YourModel) -> ComponentArray`

Default values for the model-specific exchange parameters.

```julia
defaultparams(::YourModel) = ComponentArray(; kex=1000.0, pB=0.05)
```

#### `exchangematrix(::YourModel, params, expt) -> Matrix{Float64}`

The ``N \times N`` kinetic exchange matrix. Element ``K_{ij}`` is the rate from state ``j`` to state ``i``. Column sums must be zero (conservation).

```julia
function exchangematrix(::YourModel, params, expt)
    kex = params.model.kex
    pB = params.model.pB
    pA = 1 - pB
    return [-kex*pB  kex*pA;
             kex*pB -kex*pA]
end
```

#### `populations(::YourModel, params, expt) -> Vector{Float64}`

Equilibrium populations of each state. Must sum to 1.

```julia
function populations(::YourModel, params, expt)
    pB = params.model.pB
    return [1 - pB, pB]
end
```

!!! note
    For binding models, `populations` typically depends on sample concentrations accessed via `sampleconcentrations(expt)` and the model's `moleculemap`. See `TwoStateBindingModel` for an example.

### Step 4: Include the file

Add an include statement to `src/exchange1d/models.jl`:

```julia
include("model-yourmodel.jl")
```

### Step 5: Verify

The model will automatically appear in the selection menu because `_available_models()` uses `subtypes(AbstractModel)` to discover all concrete subtypes at runtime. No further registration is needed.

### Complete Example

Here is the full implementation of the two-state exchange model (`model-twostate.jl`):

```julia
struct TwoStateModel <: AbstractModel
end

modelname(::TwoStateModel) = "Two-state exchange"
nstates(::TwoStateModel) = 2
states(::TwoStateModel) = ["A", "B"]
nmolecules(::TwoStateModel) = 1
molecules(::TwoStateModel) = Dict(:A => "observed")
defaultparams(::TwoStateModel) = ComponentArray(; kex=1000.0, pB=0.05)

function exchangematrix(::TwoStateModel, params, expt)
    kex = params.model.kex
    pB = params.model.pB
    pA = 1 - pB
    return [-kex*pB  kex*pA;
             kex*pB -kex*pA]
end

function populations(::TwoStateModel, params, expt)
    pB = params.model.pB
    return [1 - pB, pB]
end
```

## Adding a New Experiment Type

New experiment types require implementing several dispatch methods and updating the experiment classifier.

### Step 1: Create the experiment file

Create a new file `src/exchange1d/expt-yourtype.jl`.

### Step 2: Define the struct

All experiment types must include these fields:

```julia
struct YourExperiment <: AbstractExperiment
    spec::Any                                    # NMRData object
    field_teslas::Float64                         # Magnetic field strength
    sampleconcentrations::Dict{String,Float64}    # Molecule concentrations
    # Experiment-specific fields:
    your_variable::Vector{Float64}                # e.g. delays, offsets
    observed_intensities::Vector{Measurement{Float64}}
    predicted_intensities::Vector{Float64}
end
```

### Step 3: Implement required functions

#### Constructor from filename

Load and validate the NMR data, extract metadata from annotations:

```julia
function YourExperiment(filename)
    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("$filename must have annotations"))

    # Extract field strength
    field_teslas = 2π * metadata(spec, 1, :bf) /
                   gyromagneticratio(metadata(spec, 1, :nucleus))
    field_teslas = round(field_teslas; digits=2)

    # Validate and extract experiment-specific data from annotations
    your_variable = annotations(spec, :your_type, :your_key)

    observed_intensities = [0.0 ± 0.0 for _ in your_variable]
    predicted_intensities = zeros(Float64, length(your_variable))

    return YourExperiment(spec, field_teslas, sampleconcentrations(spec),
                          your_variable, observed_intensities, predicted_intensities)
end
```

#### `default_spin_params(expt::YourExperiment, nstates) -> Vector{Pair{Symbol,Any}}`

Return parameter entries needed by this experiment type. Use `field_label(expt)` to create field-specific keys:

```julia
function default_spin_params(expt::YourExperiment, nstates)
    fl = field_label(expt)
    return [:delta => fill(expt.spec[1, :offsetppm], nstates),
            Symbol("R2_", fl) => fill(10.0, nstates),
            Symbol("R1_", fl) => [1.5]]
end
```

!!! tip
    Parameters with the same key from different experiments are deduplicated — the first occurrence wins. This allows experiments at the same field to share relaxation rates automatically.

#### `default_nuisance_params(expt::YourExperiment) -> Vector{Pair{Symbol,Any}}`

Return per-experiment fitting parameters (e.g. amplitude scaling):

```julia
function default_nuisance_params(expt::YourExperiment)
    fl = field_label(expt)
    return [Symbol("YourType_", fl, "_I0") => 1.0]
end
```

#### `integrate!(expt::YourExperiment, peakppm, noiseppm, ppmwidth)`

Extract peak and noise integrals from the spectrum, updating `observed_intensities` with `Measurement` values (value ± uncertainty):

```julia
function integrate!(expt::YourExperiment, peakppm, noiseppm, ppmwidth)
    spec = expt.spec

    # Integrate noise region
    noiseselector = (noiseppm - ppmwidth / 2) .. (noiseppm + ppmwidth / 2)
    noise = std(vec(data(sum(spec[noiseselector, :]; dims=F1Dim))))

    # Integrate signal region
    signalselector = (peakppm - ppmwidth / 2) .. (peakppm + ppmwidth / 2)
    integrals = vec(data(sum(spec[signalselector, :]; dims=F1Dim)))

    # Normalise
    scale = maximum(abs, integrals)
    noise /= scale
    integrals /= scale

    return expt.observed_intensities .= integrals .± noise
end
```

#### `simulate!(expt::YourExperiment, model::AbstractModel, params::ComponentArray)`

Compute predicted intensities from the model and parameters. Update `predicted_intensities` in-place:

```julia
function simulate!(expt::YourExperiment, model::AbstractModel, params::ComponentArray)
    fl = field_label(expt)
    # Access parameters:
    #   params.model.*       — exchange parameters
    #   params.spin.*        — chemical shifts and relaxation rates
    #   params.nuisance.*    — experiment-specific parameters

    # For exchange-dependent experiments, use the Liouvillian:
    for i in eachindex(expt.your_variable)
        L = liouvillian(model, params, expt, offset, rf_field)
        # ... propagate and extract observable ...
        expt.predicted_intensities[i] = result
    end

    return nothing
end
```

#### `plot_result(expt::YourExperiment, fit_result; kwargs...)`

Create a diagnostic plot showing observed data, fitted curve, and residuals:

```julia
function plot_result(expt::YourExperiment, fit_result; kwargs...)
    x = expt.your_variable
    yobs = expt.observed_intensities
    ypred = expt.predicted_intensities

    # Upper panel: data + fit
    p1 = scatter(x, yobs; ylabel="Intensity", kwargs...)
    plot!(p1, x, ypred)

    # Lower panel: weighted residuals
    wres = (Measurements.value.(yobs) .- ypred) ./ Measurements.uncertainty.(yobs)
    p2 = scatter(x, wres; xlabel="Variable", ylabel="Residual / σ", kwargs...)

    return plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]))
end
```

### Step 4: Include the file

Add an include statement to `src/exchange1d/experiments.jl`, before the `load_experiment` function:

```julia
include("expt-yourtype.jl")
```

### Step 5: Update the experiment classifier

Add a dispatch case to `load_experiment()` in `src/exchange1d/experiments.jl`:

```julia
function load_experiment(filename)
    # ... existing code ...
    elseif "yourtype" in types
        return YourExperiment(filename)
    # ... rest of function ...
end
```

### Step 6: Update analysis registration (if needed)

If your new experiment type should be included in the automatic `analyse()` dispatch, update the matcher function in `Exchange1D.__init__()` (in `src/exchange1d/Exchange1D.jl`):

```julia
function __init__()
    rule = MultiFileRule(expts -> begin
                             oneD = filter(e -> "1d" in e.types, expts)
                             cest = filter(e -> "cest" in e.types, oneD)
                             r1cal = filter(e -> "relaxation" in e.types &&
                                                "R1" in e.features, oneD)
                             yourtype = filter(e -> "yourtype" in e.types, oneD)
                             combined = vcat(cest, r1cal, yourtype)
                             length(cest) > 0 ? combined : nothing
                         end,
                         expts -> exchange1d([e.filename for e in expts]),
                         "Exchange analysis (CEST + R1)")
    return register_analysis!(rule)
end
```

## Registering a New Analysis with `analyse()`

The `analyse()` function uses a registry-based dispatch system (see [Analysis Rules](@ref)). Exchange1D registers itself in its `__init__()` function using a `MultiFileRule`:

```julia
register_analysis!(MultiFileRule(matcher, handler, name))
```

- **matcher**: A function that receives all classified `ExperimentFile` objects and returns the subset relevant to this analysis, or `nothing` if the analysis doesn't apply
- **handler**: A function that receives the matched `ExperimentFile` objects and runs the analysis
- **name**: A label shown in the interactive menu

If you are creating an entirely new analysis module (not just extending Exchange1D), you can register it using the same pattern. See [Analysis Rules](@ref) for the full API including `SingleFileRule` for single-file analyses.
