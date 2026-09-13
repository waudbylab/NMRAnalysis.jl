# TRACT: a TROSY / anti-TROSY pair, yielding the rotational correlation time τc.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    tract()
    tract(trosy, antitrosy; tau=nothing, regions=nothing, integration=nothing,
          prompt=isinteractive())

Analyse a TRACT pair, deriving τc from the TROSY / anti-TROSY relaxation-rate difference,
then open the analysis window. The two spectra are combined into one dataset tagged by
`which ∈ {:trosy, :anti}`, sharing a single integration region.

Called with no arguments, you are asked for the two experiment folders. Relaxation delays
are taken from `tau` if given, else from each spectrum's own `vdlist`, and otherwise you
are asked for them.

# Arguments
- `tau`: relaxation delays, in seconds, one per spectrum, used for both experiments.
- `regions`: integration regions to start from, instead of the 7.5-9.5 ppm amide window.
- `integration`: a `(; peakppm, noiseppm, ppmwidth)` triple, which skips the window and
  analyses that region directly.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.
"""
function tract(trosy, antitrosy; tau=nothing, regions=nothing, integration=nothing,
               prompt::Bool=isinteractive())
    # the arguments as written, for the reproduce line
    giventrosy, givenanti = trosy, antitrosy
    trosy, antitrosy = loadspec(trosy), loadspec(antitrosy)
    ttau = @something(tau, acqusvalue(trosy, :vdlist),
                      askvector("TROSY relaxation delays", nplanesfromspec(trosy);
                                unit="s", prompt))
    atau = @something(tau, acqusvalue(antitrosy, :vdlist),
                      askvector("anti-TROSY relaxation delays",
                                nplanesfromspec(antitrosy); unit="s", prompt))

    traces = vcat(tracesfromspec(trosy), tracesfromspec(antitrosy))
    vars = vcat([(; time=Float64(t), which=:trosy) for t in ttau],
                [(; time=Float64(t), which=:anti) for t in atau])

    B0 = 2π * acqus(trosy, :bf1) / GAMMA_H
    ωN = 2π * acqus(trosy, :bf3)
    f = tractf(; B0)

    # Per-plane sources, so `series.csv` names the spectrum each row came from rather
    # than naming the TROSY one for both.
    src = vcat(fill(speclabel(trosy), length(ttau)),
               fill(speclabel(antitrosy), length(atau)))
    ds = Dataset1D(Planes(traces, vars), defaultnoisecentre(trosy),
                   speclabel(trosy), src)
    expt = isnothing(regions) ? TractExperiment(ds; ωN, f) :
           TractExperiment(ds; ωN, f, regions)
    # `tau` applies to both experiments, so it is only worth recording when the two lists
    # agree; where they differ, each spectrum's own vdlist is what reproduces the analysis.
    return run1d(expt; integration,
                 call=analysiscall("tract", giventrosy, givenanti;
                                   tau=(ttau == atau ? ttau : nothing)))
end

function tract(; kwargs...)
    println("Current directory: $(pwd())")
    trosy = askpath("TROSY experiment")
    antitrosy = askpath("anti-TROSY experiment")
    return tract(trosy, antitrosy; kwargs...)
end

# ---- 2. type ------------------------------------------------------------------

"""
    TractExperiment(dataset; ωN, f, regions=…)

Fit TROSY and anti-TROSY decays (grouped by `vars.which ∈ {:trosy, :anti}`) and derive
the rotational correlation time τc from the cross-correlated relaxation-rate difference
ΔR = R(anti) − R(trosy). `ωN` is the ¹⁵N Larmor frequency (rad s⁻¹) and `f` the
dipole/CSA cross-correlation prefactor (see `tractf`).
"""
struct TractExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    ωN::Float64
    f::Float64
end

function TractExperiment(dataset::Dataset1D; ωN, f, regions=[defaultamideregion()])
    return TractExperiment(dataset, collect(Region, regions), Float64(ωN), Float64(f))
end

"""
    defaultamideregion(; label="amide") -> Region

TRACT integrates the bulk backbone amide envelope rather than a single resolved peak, so
its default region is the conventional 7.5-9.5 ppm amide window - unlike every other
experiment's `defaultregion`, which centres a narrow window on the tallest peak.
"""
defaultamideregion(; label="amide") = Region(label, 7.5, 9.5)

# ---- 3. interface -------------------------------------------------------------

seriesmodel(::TractExperiment) = ExponentialModel()
fitaxis(::TractExperiment) = :time
groupcols(::TractExperiment) = (:which,)
primaryparam(::TractExperiment) = :tauc

# ---- 4. science ---------------------------------------------------------------

"¹H gyromagnetic ratio / rad s⁻¹ T⁻¹."
const GAMMA_H = 2.6752218744e8

# The two decay rates reach the region as :R_trosy and :R_anti (see `seriesname`), which
# is what lets η and τc sit beside them on the one row the region reports.
function postfit!(r::RegionResult, e::TractExperiment)
    haskey(r.parameters, :R_trosy) && haskey(r.parameters, :R_anti) || return nothing
    ηxy = (param(r, :R_anti) - param(r, :R_trosy)) / 2
    setpost!(r, :etaxy, ηxy)
    setpost!(r, :tauc, tracttauc(e.f, e.ωN, ηxy))
    return nothing
end

"""
    tractf(; B0, θ=17π/180) -> Float64

Dipole/CSA cross-correlation prefactor used in the TRACT τc relation, given the static
field `B0` (T). Constants follow the standard ¹⁵N–¹H amide treatment.
"""
function tractf(; B0, θ=17 * π / 180)
    μ0 = 4π * 1e-7
    γH = GAMMA_H
    γN = -2.7126180e7
    ħ = 6.62607015e-34 / 2π
    rNH = 1.02e-10
    ΔδN = 160e-6
    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = B0 * γN * ΔδN / (3 * sqrt(2))
    return p * c * (3cos(θ)^2 - 1)
end

"""
    tracttauc(f, ωN, ηxy) -> Float64

Rotational correlation time τc (ns) from the cross-correlated cross-relaxation rate
`ηxy`, by the analytic inversion of

    ηxy = f · [4·J(0) + 3·J(ωN)],   J(ω) = (2/5)·τc / (1 + ω²τc²)

i.e. `ηxy = f·(8/5·τc + 6/5·τc/(1+(ωN·τc)²))`, as in the routine this replaces.

The `(2/5)` spectral-density convention is the easy thing to lose here: an earlier version
of this docstring quoted the relation as `4/5·τc + 3/5·τc/(1+(ωN·τc)²)`, a factor of two
out from what the inversion below actually solves (the code was right, the docstring was
not). `test/analysis1d_test.jl` now round-trips this against the forward relation above, so
the two cannot drift apart again silently.
"""
function tracttauc(f, ωN, ηxy)
    x = sqrt(21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy^2 * ωN^8 + 625 * f^2 * ηxy^4 * ωN^10)
    y = cbrt(1800 * f^2 * ηxy * ωN^4 + 125 * ηxy^3 * ωN^6 + 24 * sqrt(3) * x)
    τc = (5 * ηxy) / (24 * f) -
         (336 * f^2 * ωN^2 - 25 * ηxy^2 * ωN^4) / (24 * f * ωN^2 * y) + y / (24 * f * ωN^2)
    return 1e9 * τc
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::TractExperiment) = "TRACT"

resultxfactor(e::TractExperiment) = timescale(column(dataset(e).planes, :time))[1]

function resultlabels(e::TractExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end

seriesnames(::TractExperiment) = ["TROSY", "anti-TROSY"]

function spectruminfo(::TractExperiment, vars::NamedTuple)
    which = vars.which == :trosy ? "TROSY" : "anti-TROSY"
    return "$(round(vars.time; digits=3)) s delay ($which)"
end

# :R means a relaxation rate for both TROSY and anti-TROSY decays; :etaxy and :tauc appear
# only here. The keys are ASCII even though the quantities are written η and τc, a key
# becoming a CSV column header; the typeset form is in the label.
const TRACT_PARAM_LABELS = Dict(:R => "Relaxation rate",
                                :etaxy => "CCR rate (η)",
                                :tauc => "Correlation time (τc)")
const TRACT_PARAM_UNITS = Dict(:etaxy => "s-1",
                               :tauc => "ns")

function paramlabel(::TractExperiment, name::Symbol)
    return get(TRACT_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::TractExperiment, name::Symbol)
    return get(TRACT_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
