# Nutation calibration: pulse-length calibration from a damped sinusoid.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    calibration1d(spec; durations=nothing, phase=nothing, regions=nothing,
                  integration=nothing, prompt=isinteractive())

Analyse a 1D nutation calibration, reporting the nutation frequency, 90° pulse length and
B₁ inhomogeneity, then open the analysis window.

Pulse durations are taken from `durations` if given, else from the `calibration.duration`
annotation, and otherwise you are asked for them; the modulation likewise from `phase`,
else the `calibration.model` annotation, else a question. Neither annotation is required.

# Arguments
- `durations`: pulse durations, in seconds, one per spectrum.
- `phase`: `:sine` (starting from equilibrium) or `:cosine` (from transverse
  magnetisation).
- `regions`: integration regions to start from, instead of one on the tallest peak.
- `integration`: a `(; peakppm, noiseppm, ppmwidth)` triple, which skips the window and
  analyses that region directly.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.
"""
function calibration1d(spec; durations=nothing, phase=nothing, regions=nothing,
                       integration=nothing, prompt::Bool=isinteractive())
    spec = loadspec(spec)
    t = @something(durations,
                   annotation(spec, :calibration, :duration),
                   askdurations(nplanesfromspec(spec); prompt))
    phase = @something(phase,
                       nutationphase(annotation(spec, :calibration, :model)),
                       asknutationphase(; prompt))
    ds = datasetfromspec(spec, [(; duration=Float64(d)) for d in t])
    expt = isnothing(regions) ? NutationExperiment(ds; phase) :
           NutationExperiment(ds; phase, regions)
    return run1d(expt; integration)
end

"""Pulse durations, asked for in µs (as they are set on the spectrometer) but returned in
the seconds the analysis works in."""
askdurations(n::Integer; prompt::Bool=true) = 1e-6 .* askvector("pulse durations", n;
                                                                unit="µs", prompt)

"""
    nutationphase(annotation) -> Symbol or nothing

The modulation named by a `calibration.model` annotation, or `nothing` where there is no
annotation to read - so the entry point can fall through to asking.
"""
nutationphase(::Nothing) = nothing
nutationphase(s) = string(s) == "cosine_modulated" ? :cosine : :sine

"""Offer the choice of modulation. Without prompting, a sine modulation is assumed."""
function asknutationphase(; prompt::Bool=true)
    choice = askchoice("How is the signal modulated by the nutation pulse?",
                       ["Sine:   I(t) = A sin(2πνt) exp(-Rt), from equilibrium",
                        "Cosine: I(t) = A cos(2πνt) exp(-Rt), from transverse magnetisation"];
                       prompt)
    return choice == 2 ? :cosine : :sine
end

# ---- 2. type ------------------------------------------------------------------

"""
    NutationExperiment(dataset; phase=:sine, regions=…)

Fit integrated intensity vs pulse duration (`vars.duration`) to a damped sinusoid and
derive the 90° pulse length (`1/(4ν)`) and B₁ inhomogeneity (`R/2πν`). A height
(zero-width region at the observed signal) is the natural reduction but a finite region
works identically.
"""
struct NutationExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    model::SeriesModel
end

function NutationExperiment(dataset::Dataset1D; phase::Symbol=:sine,
                            regions=[defaultregion(dataset)])
    return NutationExperiment(dataset, collect(Region, regions),
                              DampedSinusoidModel(; phase))
end

# ---- 3. interface -------------------------------------------------------------

fitaxis(::NutationExperiment) = :duration
primaryparam(::NutationExperiment) = :pulse90

# ---- 4. science ---------------------------------------------------------------

"""
    DampedSinusoidModel(; phase = :sine)

Damped sinusoid for nutation calibration: `A·sin(2π·ν·t)·exp(−R·t)` (or `cos` when
`phase = :cosine`), parameters `[A, ν, R]`. `ν` is the nutation frequency (Hz) from
which the 90° pulse length follows as `1/(4ν)`.
"""
function DampedSinusoidModel(; phase::Symbol=:sine)
    trig = phase === :cosine ? cos : sin
    function est(x, y)
        A = maximum(abs.(y))
        # rough frequency: assume ~half a period across the sampled range
        ν = 0.5 / (maximum(x) - minimum(x))
        return [A, ν, 1.0 / maximum(x)]
    end
    return CurveFitModel((x, p) -> @.(p[1] * trig(2π * p[2] * x) * exp(-p[3] * x)),
                         ["A", "ν", "R"],
                         est;
                         xlabel="Pulse duration / s")
end

# Stored in the units `PARAM_UNITS` names for them: a 90° pulse reads naturally in µs,
# never in seconds.
function postfit!(r::RegionResult, ::NutationExperiment)
    ν = param(r, :ν)
    setpost!(r, :pulse90, 1e6 / (4ν))
    setpost!(r, :inhomogeneity, 100 * param(r, :R) / (2π * ν))
    return nothing
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::NutationExperiment) = "Nutation calibration"

resultxfactor(e::NutationExperiment) = timescale(column(dataset(e).planes, :duration))[1]

function resultlabels(e::NutationExperiment)
    _, unit = timescale(column(dataset(e).planes, :duration))
    return ("Pulse duration / $unit", "Integrated intensity (a.u.)")
end

function spectruminfo(::NutationExperiment, vars::NamedTuple)
    return "$(round(1e6 * vars.duration; digits=1)) µs pulse"
end

# Own display names and units, not the shared PARAM_LABELS/PARAM_UNITS tables -
# everything about this experiment's presentation lives here. :ν, :pulse90 and
# :inhomogeneity only ever appear in this file (DampedSinusoidModel and postfit!, above).
# :R is deliberately *not* overridden here: it's this model's decay rate, not a
# relaxation rate (see the shared table's own note on why it stays bare by default).
const NUTATION_PARAM_LABELS = Dict(:ν => "Nutation frequency",
                                   :pulse90 => "90°",
                                   :inhomogeneity => "B₁ inhom.")
const NUTATION_PARAM_UNITS = Dict(:ν => " Hz",
                                  :pulse90 => " µs",
                                  :inhomogeneity => " %")

function paramlabel(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
