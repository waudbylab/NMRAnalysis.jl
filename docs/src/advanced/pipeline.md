# How an analysis works

Every analysis in NMRAnalysis.jl, 1D or 2D, is the same three stages. Most analyses use
only the first two, and some only the first.

## Stage 1: measure

The user picks the things to be measured: regions in 1D, peaks in 2D. Each one is reduced
to a single number in every plane. In 1D that means integrating the region; in 2D it means
fitting a lineshape, with overlapping peaks fitted together as a cluster, which changes how
the fit is run but not what comes out of it.

The result is a **series** for each region or peak: one measured value per plane, together
with whatever else the measurement found, which in 2D is the peak's position and linewidth
in each plane. The experiment supplies what distinguishes one plane from another, the
relaxation delay or the gradient strength or the saturation offset, and labels the planes
accordingly.

Where the same region or peak is measured under more than one condition, each condition is
its own series: TROSY and anti-TROSY in TRACT, reference and saturated in a heteronuclear
NOE, isotropic and aligned in an RDC measurement, and one series per experiment in a
combined exchange analysis.

## Stage 2: fit

Each region or peak's series are fitted to give the numbers actually reported for it: a
relaxation rate, a diffusion coefficient and hydrodynamic radius, a 90° pulse length, an
NOE ratio, τc from a TROSY/anti-TROSY pair, a scalar coupling and an RDC from the four
planes of an aligned pair.

This stage sees all of an entity's series at once, so anything that combines conditions
belongs here rather than being a special case somewhere else.

## Stage 3: global fit

The series and the per-entity results together are fitted for parameters shared by every
entity: a dissociation constant across a titration, an exchange rate and populations across
a set of exchange experiments.

A global fit may also write back per-entity numbers, because determining the shared
parameter is what determines them. A titration's free and bound shifts, and the chemical
shift perturbation derived from them, come out of the global fit but describe one peak.

## What goes where

Each stage has an output file, and a number is written to the file matching what it
describes rather than the stage that produced it.

| scope | file | examples |
|---|---|---|
| varies from plane to plane | `series.csv` | integrals, amplitudes, fitted positions and linewidths |
| describes one region or peak | `results.csv` | relaxation rates, τc, D and r_H, NOE ratios, CSPs |
| shared by every region or peak | `global.csv` | a titration K_d, an exchange k_ex |

Kinetics stops after stage 1: the intensity of each region against time is the answer, and
there is no `results.csv`. A relaxation or diffusion measurement stops after stage 2 and
has no `global.csv`. A combined exchange analysis has a single region and nothing to fit
per region, so its stage 2 is empty and everything is global.

See [Output and Interface Conventions](conventions.md) for the layout and column rules of
these files.

## Where each stage lives in the code

| stage | `Analysis1D` | `GUI2D` |
|---|---|---|
| 1, measure | `integrate(region, expt)` | `fit!(cluster, expt)` |
| 2, fit | `postfit!(results, expt)` | `postfit!(peak, expt)` |
| 3, global fit | `postfitglobal!(results, expt)` | `postfitglobal!(expt)` |

Stage 2 receives the whole region or peak, not one series or plane at a time, which is
what lets a quantity combining conditions (TRACT's τc) be an ordinary post-fit. In
`Analysis1D` a region's fitted parameters are named apart by the series they came from
(`:R_trosy`, `:R_anti`) and kept in one set with everything derived from them, so one
region is one row of results.
