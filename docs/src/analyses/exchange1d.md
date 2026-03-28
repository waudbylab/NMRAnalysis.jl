# 1D Exchange Analysis

The Exchange1D module in NMRAnalysis.jl provides tools for analysing chemical exchange processes using 1D NMR experiments. It performs joint fitting of CEST (Chemical Exchange Saturation Transfer) and R1 relaxation data using Bloch-McConnell equations to extract exchange kinetics, populations, and chemical shift differences.

## Launching Exchange Analysis

### Via automatic dispatch

If your experiments have appropriate annotations, the `analyse()` function will detect CEST and R1 experiments and offer exchange analysis automatically:

```julia
using NMRAnalysis
analyse(["data/101", "data/102", "data/103"])
```

See [Automatic Analysis Dispatch](@ref) below for details on how `analyse()` works.

### Direct launch with a directory

Call `exchange1d()` with a folder path to select experiments interactively:

```julia
exchange1d("data/exchange_expts")
```

### Direct launch with specific files

Provide experiment paths directly:

```julia
exchange1d(["data/101", "data/102", "data/103", "data/104"])
```

### Direct launch with a file dialog

Call with no arguments to open a file selection dialog:

```julia
exchange1d()
```

## Supported Experiment Types

Exchange1D currently supports two experiment types, which can be combined in a joint fit:

| Experiment | Annotation type | Annotation feature | Description |
|---|---|---|---|
| **CEST** | `"cest"` | — | Saturation transfer profiles at varying offsets |
| **R1 relaxation** | `"relaxation"` | `"R1"` | Longitudinal relaxation decay or inversion-recovery |

At least one CEST experiment is required. R1 experiments are optional and help constrain the longitudinal relaxation rate during fitting.

## Available Models

Four exchange models are available, selected interactively at the start of the analysis:

| Model | States | Parameters | Use case |
|---|---|---|---|
| **No exchange** | 1 | — | Null model (no exchange) |
| **Two-state exchange** | 2 | `kex`, `pB` | Intramolecular exchange between two conformations |
| **Three-state exchange** | 3 | `koffB`, `pB`, `koffC`, `pC` | Linear three-state exchange (A ⇌ B, A ⇌ C) |
| **Two-state binding** | 2 | `Kd`, `koff` | Bimolecular binding; requires sample concentrations |

!!! tip
    The two-state binding model uses the quadratic binding equation to calculate populations from the dissociation constant (`Kd`) and total concentrations of the observed molecule and titrant. These concentrations are read from the sample metadata in the NMR annotations.

## Analysis Workflow

The interactive text-based workflow proceeds through the following steps:

### 1. Model selection

A menu is presented listing all available exchange models. Select the model appropriate for your system.

### 2. Experiment loading

Each input file is loaded and classified automatically based on its annotations. The loaded experiments are listed with their type and magnetic field strength.

### 3. Molecule mapping (binding models only)

For models involving multiple molecules (e.g. two-state binding), you are prompted to assign molecule roles (observed species, titrant) to the sample components found in the experiment metadata.

### 4. Peak integration

You are prompted for three values:
- **Peak position** (ppm): Centre of the signal of interest
- **Integration width** (ppm): Width of the integration region (default: 0.1 ppm)
- **Noise position** (ppm): Centre of a signal-free region for noise estimation

All experiments are integrated at these positions. Noise is estimated from the standard deviation of the noise region across the variable dimension (e.g. saturation offsets for CEST, delays for R1).

### 5. Parameter review and fitting

Default parameters are assembled from the model and experiments. An interactive editor lets you review and modify initial values before fitting. After fitting, results are displayed as a table and diagnostic plots.

You can then choose to:
- **Save** the results to a folder
- **Adjust** parameters and refit
- **Quit** without saving

### 6. Saving results

If you choose to save, you are prompted for an output folder. The following files are written:
- `exchange1d_fit.pdf` — Combined plot of all experiments
- `exchange1d_expt_N.pdf` — Individual experiment plots
- `exchange1d_params.txt` — Fitted parameters with uncertainties

## Parameter Structure

Parameters are organised into three sections using a `ComponentArray`:

| Section | Contents | Example keys |
|---|---|---|
| **model** | Exchange kinetics and populations | `kex`, `pB`, `Kd`, `koff` |
| **spin** | Chemical shifts and field-dependent relaxation rates | `delta`, `R2_14p1T`, `R1_14p1T` |
| **nuisance** | Per-experiment amplitude and correction factors | `R1_14p1T_I0`, `R1_14p1T_inv_factor` |

Field-dependent parameters are tagged with the magnetic field strength (e.g. `R2_14p1T` for 14.1 T). When experiments at multiple fields are combined, separate parameters are created for each field.

Chemical shifts (`delta`) and transverse relaxation rates (`R2`) are vectors with one value per state. Longitudinal relaxation (`R1`) is shared across states (length 1).

## Return Value

The `exchange1d()` function returns a `FitResult` with the following fields:

| Field | Type | Description |
|---|---|---|
| `params` | `ComponentArray{Measurement}` | Fitted parameters with uncertainties |
| `params_value` | `ComponentArray{Float64}` | Fitted parameter values (no uncertainties) |
| `params0` | `ComponentArray{Float64}` | Initial parameters used for the fit |
| `chi2` | `Float64` | Chi-squared statistic |
| `reduced_chi2` | `Float64` | Reduced chi-squared (``\chi^2 / \text{dof}``) |
| `cov` | `Matrix{Float64}` | Parameter covariance matrix |
| `nobs` | `Int` | Number of data points |
| `nparams` | `Int` | Number of fitted parameters |
| `dof` | `Int` | Degrees of freedom |
| `prob` | `ExchangeProblem` | Reference to the fitted problem |

## Automatic Analysis Dispatch

The `analyse()` function provides a convenient way to automatically detect and run appropriate analysis routines based on experiment metadata.

### How it works

1. **Classification**: Each input file is loaded and classified by its `types` (e.g. `"1d"`, `"cest"`, `"relaxation"`) and `features` (e.g. `"R1"`, `"nutation"`)
2. **Matching**: Registered analysis rules are matched against the classified files
3. **Selection**: If multiple analyses match, an interactive menu is presented
4. **Execution**: Selected analyses are run and results returned

### Usage

```julia
# Single file
result = analyse("data/101")

# Multiple files
results = analyse(["data/101", "data/102", "data/103"])
```

When exchange analysis is available (at least one CEST experiment is present), it appears in the menu as "Exchange analysis (CEST + R1)".

!!! note
    The `analyse()` function is extensible. Analysis modules register themselves at load time, so all available analyses are automatically discovered. See [Analysis Rules](@ref) for details on how to register custom analysis routines.

## Theoretical Background

### Bloch-McConnell Equations

Exchange analysis uses the Bloch-McConnell formalism, which extends the Bloch equations to multiple exchanging states. For ``N`` states, the magnetisation vector ``\mathbf{M}`` evolves as:

```math
\frac{d\mathbf{M}}{dt} = \mathbf{L} \, \mathbf{M}
```

where ``\mathbf{L}`` is the ``3N \times 3N`` Liouvillian superoperator incorporating:
- Chemical shift evolution (``\Omega_i`` for each state)
- Relaxation (``R_1``, ``R_{2,i}`` for each state)
- RF fields (spin-lock or saturation)
- Chemical exchange (kinetic rate matrix ``\mathbf{K}``)

### CEST Simulation

For CEST experiments, the Liouvillian is augmented with an inhomogeneous term to account for relaxation back to equilibrium. The predicted CEST profile is computed by propagating the initial equilibrium magnetisation through the saturation period:

```math
\mathbf{M}(T_\text{sat}) = \exp(\mathbf{L}_\text{inhom} \cdot T_\text{sat}) \, \mathbf{M}_0
```

The observed intensity is the sum of ``M_z`` components across all states.

### R1 Simulation

R1 experiments are fitted independently of the exchange model (R1 decay is not sensitive to chemical exchange under typical conditions). Depending on the experiment type:

- **Exponential decay**: ``I(t) = I_0 \exp(-R_1 t)``
- **Inversion recovery**: ``I(t) = I_0 (1 - f \exp(-R_1 t))`` where ``f`` is the inversion factor

### Exchange Matrix

The kinetic exchange matrix ``\mathbf{K}`` encodes the rates of interconversion between states. Each element ``K_{ij}`` is the rate from state ``j`` to state ``i``, and column sums are zero (conservation of total magnetisation).

For the two-state model:

```math
\mathbf{K} = \begin{pmatrix} -k_\text{ex} p_B & k_\text{ex} p_A \\ k_\text{ex} p_B & -k_\text{ex} p_A \end{pmatrix}
```
