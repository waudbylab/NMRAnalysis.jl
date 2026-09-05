# 1D Exchange Analysis

The Exchange1D module in NMRAnalysis.jl provides tools for analysing chemical exchange processes using 1D NMR experiments. It performs joint fitting of CEST (Chemical Exchange Saturation Transfer), on- and off-resonance R1ρ relaxation dispersion, and R1 relaxation data using Bloch-McConnell equations to extract exchange kinetics, populations, and chemical shift differences.

This page covers how to launch an analysis and what to expect from the interactive workflow. For the exchange models themselves — what each one represents physically and which parameters it exposes — see [Exchange Models](models.md). For the underlying Bloch-McConnell and simulation equations, see [Theory and Calculation Methods](theory.md).

!!! note "Relationship to `r1rho()`"
    This is a different tool from the standalone [`r1rho()`](../r1rho.md) GUI. `r1rho()` fits
    on-resonance R1ρ dispersion data alone, with its own dedicated interactive interface.
    `exchange1d()` instead performs joint Bloch-McConnell fitting across any combination of
    CEST, R1, and on-/off-resonance R1ρ experiments against a shared exchange model, through
    a text-based workflow rather than a GUI. Use `r1rho()` for routine on-resonance dispersion
    measurements; use `exchange1d()` when you need to combine multiple experiment types, or
    both resonance offsets, in a single global fit.

## Launching Exchange Analysis

### Via automatic dispatch

If your experiments have appropriate annotations, the `analyse()` function will detect a CEST or off-resonance R1ρ experiment (any accompanying R1 or on-resonance R1ρ experiments are included automatically) and offer exchange analysis:

```julia
using NMRAnalysis
analyse(["data/101", "data/102", "data/103"])
```

See [Automatic Analysis](../analyse.md) for details on how dispatch works.

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

Exchange1D currently supports four experiment types, which can be combined in a joint fit:

| Experiment | Annotation type | Annotation feature | Description |
|---|---|---|---|
| **CEST** | `"cest"` | — | Saturation transfer profiles at varying offsets |
| **R1ρ (on-resonance)** | `"r1rho"` | `"on_resonance"` | Relaxation dispersion vs spin-lock field strength, spin-lock applied on-resonance |
| **R1ρ (off-resonance)** | `"r1rho"` | `"off_resonance"` | Relaxation dispersion vs saturation offset, at a fixed spin-lock field strength |
| **R1 relaxation** | `"relaxation"` | `"R1"` | Longitudinal relaxation decay or inversion-recovery |

When using automatic dispatch via `analyse()`, at least one CEST or off-resonance R1ρ
experiment must be present to offer exchange analysis; on-resonance R1ρ and R1 experiments
can be included in the joint fit but do not by themselves trigger it. Calling `exchange1d()`
directly accepts any combination of the four experiment types above; R1 experiments are
always optional and help constrain the longitudinal relaxation rate during fitting.

## Analysis Workflow

The interactive text-based workflow proceeds through the following steps:

### 1. Model selection

A menu is presented listing all available exchange models. Select the model appropriate for
your system — see [Exchange Models](models.md) for what each one represents and when to use it.

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
| **model** | Exchange kinetics and populations (see [Exchange Models](models.md)) | `logkex`, `dGB`, `logKd`, `logkoff` |
| **spin** | Chemical shifts and field-dependent relaxation rates | `delta`, `R2_14p1T`, `R1_14p1T` |
| **nuisance** | Per-experiment amplitude and correction factors | `R1_14p1T_I0`, `R1_14p1T_inv_factor` |

Field-dependent parameters are tagged with the magnetic field strength (e.g. `R2_14p1T` for 14.1 T). When experiments at multiple fields are combined, separate parameters are created for each field.

Chemical shifts (`delta`) and transverse relaxation rates (`R2`) are vectors with one value per state. Longitudinal relaxation (`R1`) is shared across states (length 1).

!!! note "Scripting against the fit"
    `exchange1d()` also returns a `FitResult` object, for scripted analyses that need
    programmatic access to the fitted parameters, uncertainties, and fit statistics rather
    than (or in addition to) the saved files above. See the `FitResult` entry in the
    [Exchange1D API reference](../../api/exchange1d.md) for its fields.

## Automatic dispatch

When exchange analysis is available (at least one CEST or off-resonance R1ρ experiment is present), it appears in the `analyse()` menu as "Exchange analysis (CEST / R1rho)". See [Automatic Analysis](../analyse.md) for how dispatch works in general, and [Analysis Rules](../../advanced/analysis_rules.md) for how to register custom analysis routines.
