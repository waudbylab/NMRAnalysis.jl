# Nutation calibration: pulse-length calibration from a damped sinusoid.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    calibration1d(spec; durations=nothing, phase=nothing, regions=nothing, integration=nothing)

Analyse a 1D nutation calibration, reporting the nutation frequency, 90° pulse length and
B₁ inhomogeneity. Pulse durations come from `durations`, else the `calibration.duration`
annotation. The modulation defaults to the `calibration.model` annotation
(`"cosine_modulated"` selects a cosine); pass `phase=:sine`/`:cosine` to override.
"""
function calibration1d(spec; durations=nothing, phase=nothing, regions=nothing,
                       integration=nothing)
    spec = loadspec(spec)
    t = something(durations, annotation(spec, :calibration, :duration))
    ph = isnothing(phase) ?
         (annotation(spec, :calibration, :model) == "cosine_modulated" ? :cosine : :sine) : phase
    ds = datasetfromspec(spec, [(; duration=Float64(d)) for d in t])
    expt = isnothing(regions) ? NutationExperiment(ds; phase=ph) :
           NutationExperiment(ds; phase=ph, regions)
    return run1d(expt; integration)
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
    setpost!(r, :inhomogeneity, param(r, :R) / (2π * ν))
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
                                  :pulse90 => " µs")

function paramlabel(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
