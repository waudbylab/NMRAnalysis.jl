# Nutation calibration: pulse-length calibration from a damped sinusoid, over one power
# level or several.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

"""
    calibration1d(spec; durations=nothing, phase=nothing, power=nothing, regions=nothing,
                  integration=nothing, window=true, prompt=isinteractive())
    calibration1d(specs::AbstractVector; ...)

Analyse a 1D nutation calibration, reporting the nutation frequency, 90° pulse length and
B₁ inhomogeneity, then open the analysis window.

Given several spectra recorded at different power levels, each is fitted separately and the
results combined into a calibration curve: the field at a reference power and the
`linearity` of the amplifier's response, 1 being the ideal ν₁ ∝ √W. A single spectrum still
calibrates, with the linearity assumed ideal. Either way
[`B1Calibration`](@ref)`(results)` turns the analysis into the calibration that `exchange1d`
takes, so CEST and R₁ρ are simulated with measured field strengths and a measured B₁
inhomogeneity rather than nominal ones.

Pulse durations are taken from `durations` if given, else from the `calibration.duration`
annotation, and otherwise you are asked for them; the modulation likewise from `phase`,
else the `calibration.model` annotation, else a question; the power level from `power`, else
the `calibration.power` annotation, else a question. None of the annotations is required,
and an unknown power costs only the calibration curve, not the analysis.

# Arguments
- `durations`: pulse durations, in seconds, one per spectrum of the series. Applies to every
  spectrum given; where they differ, each spectrum's own annotation is what to rely on.
- `phase`: `:sine` (starting from equilibrium) or `:cosine` (from transverse
  magnetisation).
- `power`: the power each spectrum was recorded at, as `Power`s or as numbers in dB, one per
  spectrum.
- `regions`: integration regions to start from, instead of one on the tallest peak.
- `integration`: a `(; peakppm, noiseppm, ppmwidth)` triple, which skips the window and
  analyses that region directly.
- `window`: whether to open the analysis window at all. `false` analyses the default region
  (the tallest peak) without one, which is how `B1Calibration(specs)` fits a calibration
  without interrupting the analysis that asked for it.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.

# Example
```julia
results = calibration1d(["cal/1", "cal/2", "cal/3"])
cal = B1Calibration(results)      # ready for exchange1d(...; calibration=cal)
```
"""
function calibration1d(specs::AbstractVector; durations=nothing, phase=nothing, power=nothing,
                       regions=nothing, integration=nothing, window::Bool=true,
                       prompt::Bool=isinteractive())
    isempty(specs) && throw(ArgumentError("no calibration experiments given"))
    given = specs                      # the argument as written, for the reproduce line
    specs = loadspec.(specs)
    # The modulation is a property of the pulse sequence, so it is read once, from the first
    # spectrum, rather than asked for again per power level.
    phase = @something(phase,
                       nutationphase(annotation(first(specs), :calibration, :model)),
                       asknutationphase(; prompt))
    t = [@something(durations,
                    annotation(spec, :calibration, :duration),
                    askdurations(nplanesfromspec(spec); prompt)) for spec in specs]
    dB = @something(powerdb(power),
                    powerdb(annotation.(specs, :calibration, :power)),
                    askpowers(length(specs); prompt),
                    Some(nothing))
    isnothing(dB) || length(dB) == length(specs) ||
        throw(ArgumentError("got $(length(dB)) power levels for $(length(specs)) " *
                            "experiments - one power level each, in dB"))

    traces = reduce(vcat, tracesfromspec.(specs))
    vars = reduce(vcat, [planevars(t[i], isnothing(dB) ? nothing : dB[i])
                         for i in eachindex(specs)])
    src = reduce(vcat, [fill(speclabel(specs[i]), length(t[i])) for i in eachindex(specs)])
    ds = Dataset1D(Planes(traces, vars), defaultnoisecentre(first(specs)),
                   speclabel(first(specs)), src)

    expt = isnothing(regions) ? NutationExperiment(ds; phase) :
           NutationExperiment(ds; phase, regions)
    window || return analyse(expt)
    # `durations` applies to every spectrum, so it is only worth recording when the lists
    # agree; where they differ, each spectrum's own annotation reproduces the analysis.
    return run1d(expt; integration,
                 call=analysiscall("calibration1d", given;
                                   durations=(allequal(t) ? first(t) : nothing), phase,
                                   power=dB))
end

function calibration1d(spec; kwargs...)
    return calibration1d([spec]; kwargs...)
end

"""Per-plane variables for one spectrum of a calibration: the pulse duration always, and the
power level where it is known. Powers are rounded to the 0.01 dB they are set to on the
spectrometer, so that the series they name carries a legible column header."""
function planevars(durations, dB)
    isnothing(dB) && return [(; duration=Float64(d)) for d in durations]
    return [(; duration=Float64(d), power=round(Float64(dB); digits=2)) for d in durations]
end

"""Pulse durations, asked for in µs (as they are set on the spectrometer) but returned in
the seconds the analysis works in."""
function askdurations(n::Integer; prompt::Bool=true)
    return 1e-6 .* askvector("pulse durations", n;
                             unit="µs", prompt)
end

"""
    powerdb(power) -> Vector{Float64} or nothing

Power levels as dB attenuation: `Power`s converted, bare numbers taken to be dB already, and
`nothing` where any of them is missing - a calibration curve needs every power or none.
"""
powerdb(::Nothing) = nothing
powerdb(p) = powerdb([p])
function powerdb(powers::AbstractVector)
    any(isnothing, powers) && return nothing
    return [p isa Power ? db(p) : Float64(p) for p in powers]
end

"""Power levels, asked for in the dB attenuation the spectrometer sets them as. Returns
`nothing` without prompting, an unknown power costing only the calibration curve."""
function askpowers(n::Integer; prompt::Bool=true)
    prompt || return nothing
    return askvector("pulse power levels", n; unit="dB", prompt)
end

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
                       ["Sine:   I(t) = A sin(2πνt) exp(-½(2πσνt)²), from equilibrium",
                        "Cosine: I(t) = A cos(2πνt) exp(-½(2πσνt)²), from transverse magnetisation"];
                       prompt)
    return choice == 2 ? :cosine : :sine
end

# ---- 2. type ------------------------------------------------------------------

"""
    NutationExperiment(dataset; phase=:sine, regions=…)

Fit integrated intensity vs pulse duration (`vars.duration`) to a damped sinusoid and
derive the 90° pulse length (`1/(4ν)`) and B₁ inhomogeneity. Where the planes carry more
than one power level (`vars.power`, in dB), each is fitted as its own series and the results
combined into a calibration curve. A height (zero-width region at the observed signal) is
the natural reduction but a finite region works identically.
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

# One series per power level, but only where there is more than one: a single-power
# calibration reports `nu` and `pulse90`, not `nu_-12.34` and `pulse90_-12.34`.
function groupcols(e::NutationExperiment)
    p = dataset(e).planes
    hasvar(p, :power) && !allequal(column(p, :power)) && return (:power,)
    return ()
end

# ---- 4. science ---------------------------------------------------------------

"""
    DampedSinusoidModel(; phase = :sine)

Nutation of an inhomogeneous B₁ field: `A·sin(2π·ν·t)·exp(−½(2π·σ·ν·t)²)` (or `cos` when
`phase = :cosine`), parameters `[A, ν, σ]`.

`ν` is the mean nutation frequency (Hz), from which the 90° pulse length follows as
`1/(4ν)`, and `σ` is the *fractional* width of the B₁ distribution, which is what damps the
nutation: averaging `sin(2πνt)` over a Gaussian distribution of ν with mean ν̄ and standard
deviation σν̄ gives `sin(2πν̄t)·exp(−½(2πσν̄t)²)`, a Gaussian decay rather than an
exponential one. Fitting the envelope this way makes σ the same quantity
[`B1Distribution`](@ref) later samples, instead of an exponential rate needing conversion.

Relaxation during the pulse damps the nutation too, and is not separated here: σ is
therefore an upper bound, worst at the lowest power where the pulses are longest, which is
why a multi-power calibration reports the smallest of its estimates.

Only σ² enters, so the fitted sign is arbitrary and `postfit!` reports |σ|.
"""
function DampedSinusoidModel(; phase::Symbol=:sine)
    trig = phase === :cosine ? cos : sin
    est(x, y) = [maximum(abs.(y)), estimatefrequency(x, y, trig), 0.05]
    return CurveFitModel((x, p) -> @.(p[1] * trig(2π * p[2] * x) *
                                      exp(-0.5 * (2π * p[3] * p[2] * x)^2)),
                         ["A", "nu", "sigma"],
                         est;
                         xlabel="Pulse duration / s")
end

"""
    estimatefrequency(x, y, trig) -> Float64

Initial estimate of the nutation frequency: the frequency whose modulation correlates best
with the data, scanned over the frequencies the sampling can resolve.

How much of a period a calibration covers is not something the data say - the annotated
sequence runs to a nominal 720°, and the B₁ inhomogeneity is only measurable over several
periods - so assuming half a period, as this did before, starts the fit a factor of four
out on a real experiment. A coarse matched filter costs nothing next to the fit it seeds,
and the scan runs from half a period across the sampled range up to the Nyquist frequency
of the pulse-duration increment.
"""
function estimatefrequency(x, y, trig)
    span = maximum(x) - minimum(x)
    length(x) < 3 && return 0.5 / span
    nyquist = 0.5 * (length(x) - 1) / span
    νs = range(0.5 / span, nyquist; length=200)
    return argmax(ν -> abs(sum(y .* trig.(2π * ν .* x))), νs)
end

# Stored in the units `PARAM_UNITS` names for them: a 90° pulse reads naturally in µs, never
# in seconds, and an inhomogeneity as a percentage.
function postfit!(r::RegionResult, e::NutationExperiment)
    p = dataset(e).planes
    # With one power level the group is empty, so the power is not in the parameter names
    # and is recorded here instead: a 90° pulse means nothing without it.
    if hasvar(p, :power) && isempty(first(r.series).group)
        setpost!(r, :power, first(column(p, :power)))
    end
    for s in r.series
        ν = param(r, seriesname(:nu, s.group))
        setpost!(r, seriesname(:pulse90, s.group), 1e6 / (4ν))
    end

    # The smallest estimate across power levels: every one is inflated by relaxation during
    # the pulse, least so at high power where the pulses are shortest.
    σ = [abs(param(r, seriesname(:sigma, s.group))) for s in r.series]
    setpost!(r, :inhomogeneity, 100 * σ[argmin(Measurements.value.(σ))])

    # The calibration curve, where several power levels were measured. `nu1ref` is the
    # field at `powerref` and `linearity` the fitted exponent relative to the ideal power
    # law; with one power level they would only repeat `nu` and `power`.
    if hasvar(p, :power) && length(r.series) > 1
        cal = B1Calibration(r)
        setpost!(r, :powerref, db(refpower(cal)))
        setpost!(r, :nu1ref, ν1ref(cal))
        setpost!(r, :linearity, linearity(cal))
    end
    return nothing
end

"""
    B1Calibration(results::AbstractVector{RegionResult}; kwargs...)
    B1Calibration(result::RegionResult; kwargs...)
    B1Calibration(specs::AbstractVector; kwargs...)
    B1Calibration(spec; kwargs...)

The B₁ calibration a nutation analysis measured: the field fitted at each power level, and
the smallest of its B₁ inhomogeneity estimates.

Pass the results of [`calibration1d`](@ref), having checked the fits in the analysis window,
or pass the calibration experiments themselves to have them fitted without a window opening
(which is what `exchange1d(…; calibration=["cal/1", …])` does; keywords are forwarded to
`calibration1d`). Where several regions were integrated, the first is the calibration; the
others are presumably there for comparison.

Every power level has to be known, since a calibration is a map from power to field: give
them with `power=`, or annotate the sequence with `calibration.power`.

# Example
```julia
cal = B1Calibration(calibration1d(["cal/1", "cal/2"]))
hz(Power(12.0, :dB), cal)     # field at a power the calibration did not measure
```
"""
function B1Calibration(r::RegionResult; nuc=nothing, source::AbstractString="")
    isempty(r.series) && throw(ArgumentError("region \"$(r.region)\" has no fitted series"))
    haskey(r.parameters, seriesname(:nu, first(r.series).group)) ||
        throw(ArgumentError("region \"$(r.region)\" has no fitted nutation frequency; " *
                            "was the analysis run with fitting switched off?"))
    haskey(r.parameters, :power) || !isempty(first(r.series).group) ||
        throw(ArgumentError("the power level of the calibration is unknown, so it cannot " *
                            "be turned into a calibration: pass power=… to calibration1d, " *
                            "or annotate the sequence with calibration.power"))
    power = Power{Float64}[]
    ν1 = Measurement{Float64}[]
    for s in r.series
        dB = isempty(s.group) ? param(r, :power) : s.group.power
        push!(power, Power(Float64(dB), :dB))
        push!(ν1, param(r, seriesname(:nu, s.group)))
    end
    σ = [abs(param(r, seriesname(:sigma, s.group))) for s in r.series]
    return B1Calibration(power, ν1; nuc,
                         inhomogeneity=Measurements.value(σ[argmin(Measurements.value.(σ))]),
                         source)
end

function B1Calibration(results::AbstractVector{RegionResult}; kwargs...)
    isempty(results) && throw(ArgumentError("no regions were analysed"))
    return B1Calibration(first(results); kwargs...)
end

function B1Calibration(specs::AbstractVector; prompt::Bool=isinteractive(), kwargs...)
    results = calibration1d(specs; prompt, window=false, kwargs...)
    return B1Calibration(results; source=join(string.(specs), ", "))
end

function B1Calibration(spec::Union{AbstractString,Integer}; kwargs...)
    return B1Calibration([spec]; kwargs...)
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::NutationExperiment) = "Nutation calibration"

resultxfactor(e::NutationExperiment) = timescale(column(dataset(e).planes, :duration))[1]

function resultlabels(e::NutationExperiment)
    _, unit = timescale(column(dataset(e).planes, :duration))
    return ("Pulse duration / $unit", "Integrated intensity (a.u.)")
end

function spectruminfo(::NutationExperiment, vars::NamedTuple)
    info = "$(round(1e6 * vars.duration; digits=1)) µs pulse"
    haskey(vars, :power) && return "$info at $(vars.power) dB"
    return info
end

# Own display names and units, not the shared PARAM_LABELS/PARAM_UNITS tables -
# everything about this experiment's presentation lives here. The keys are ASCII (:nu, not
# :ν) because they become CSV column headers; the typeset names are in the labels.
#
# :sigma is the fitted fractional width of the B₁ distribution and :inhomogeneity the
# smallest of those estimates as a percentage - the same quantity in different units, which
# is why they are named apart rather than sharing a name with two units.
const NUTATION_PARAM_LABELS = Dict(:nu => "Nutation frequency",
                                   :pulse90 => "90°",
                                   :sigma => "B₁ inhom. (σ)",
                                   :inhomogeneity => "B₁ inhom.",
                                   :power => "Power",
                                   :powerref => "Reference power",
                                   :nu1ref => "ν₁ at reference power",
                                   :linearity => "Linearity")
const NUTATION_PARAM_UNITS = Dict(:nu => "Hz",
                                  :pulse90 => "us",
                                  :inhomogeneity => "%",
                                  :power => "dB",
                                  :powerref => "dB",
                                  :nu1ref => "Hz")

function paramlabel(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
