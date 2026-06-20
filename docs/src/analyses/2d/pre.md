# Paramagnetic Relaxation Enhancements (PREs)

!!! warning "Under development"
    PRE fitting in NMRAnalysis.jl is under active development. Results should be
    validated carefully before use in publication.

Paramagnetic relaxation enhancement (PRE) experiments measure the distance-dependent
acceleration of nuclear relaxation caused by a paramagnetic centre (e.g. a spin label,
metal ion, or paramagnetic cosolvent). The PRE rate Γ (s⁻¹ per unit concentration)
is extracted by simultaneously fitting peak lineshapes across a series of spectra recorded
at different paramagnetic agent concentrations.

The model accounts for two effects of the paramagnetic agent:

- **Linebroadening**: the direct-dimension transverse relaxation rate R2x increases by
  `Γ × conc`; in HMQC experiments the indirect-dimension rate R2y is also broadened.
- **Intensity reduction**: the peak amplitude decreases as `exp(−Γ × conc × Trelax)`.

<!-- screenshot: assets/pre2d.png -->

## Usage

### Protein PRE (diamagnetic / paramagnetic pair)

For a comparison of diamagnetic and paramagnetic protein samples (e.g. with a MTSL spin
label), set concentrations to `[0, 1]`:

```julia
using NMRAnalysis

pre2d(
    ["diamagnetic/pdata/1", "paramagnetic/pdata/1"],
    [0, 1],
    :hsqc,
    0.01
)
```

### Solvent PRE (concentration series)

For a titration of a paramagnetic cosolvent (e.g. Gd-DTPA or DSS), provide the actual
concentrations:

```julia
pre2d(
    ["0mM/pdata/1", "1mM/pdata/1", "5mM/pdata/1", "10mM/pdata/1"],
    [0.0, 1.0, 5.0, 10.0],
    :hmqc,
    0.0089
)
```

## Parameters

- **`expttype`** (`:hsqc` or `:hmqc`): In an HSQC experiment, the paramagnetic agent
  broadens only the directly detected ¹H linewidth (R2x). In an HMQC experiment, the
  heteronuclear coherence during the INEPT transfer is also broadened, so the indirect
  dimension linewidth (R2y) is additionally affected.

- **`Trelax`**: The total magnetisation transfer time during which PRE-induced relaxation
  occurs. This depends on the pulse sequence and must be set carefully. For a standard
  HSQC, it is twice the ``1/(4J)`` INEPT delay; for HMQC it includes additional coherence
  transfer periods.

## Output

Clicking **Save to folder** writes all results to `results.csv`. The derived
column is:

| Column | Description |
|--------|-------------|
| `PRE`, `PRE_err` | Fitted PRE rate Γ (s⁻¹ per unit concentration) and uncertainty |

See [Peak Lists and Output Files](peaklistformats.md) for the full format. Plot Γ
against residue number with `summaryplot(expt)`.
