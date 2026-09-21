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
  (the tallest peak) without one, which is how [`calibrationanalysis`](@ref) fits a
  calibration without interrupting the analysis that asked for it.
- `prompt`: whether to ask for anything that could not be determined. Defaults to `false`
  outside an interactive session, where a missing value raises an error instead.

# Example
```julia
results = calibration1d(["cal/1", "cal/2", "cal/3"])
cal = B1Calibration(results)      # ready for exchange1d(...; calibration=cal)
```
"""
function calibration1d(specs::AbstractVector; integration=nothing, window::Bool=true,
    kwargs...)
    expt, _, _, call = nutationexperiment(specs; kwargs...)
    window && return run1d(expt; integration, call)
    return analyse(expt)
end

function calibration1d(spec; kwargs...)
    return calibration1d([spec]; kwargs...)
end

"""
    nutationexperiment(specs; durations, phase, power, regions, prompt)
        -> (expt, dataset, regions, call)

Everything [`calibration1d`](@ref) resolves before it fits: the spectra loaded, their pulse
durations and power levels found, the experiment built from them, and the call that would
repeat the analysis.

Separate from the entry point so that a calibration fitted on behalf of another analysis
can be saved afterwards from the same pieces the entry point would have saved - see
[`calibrationanalysis`](@ref).
"""
function nutationexperiment(specs::AbstractVector; durations=nothing, phase=nothing,
    power=nothing, regions=nothing,
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
    vars = reduce(vcat,
        [planevars(t[i], isnothing(dB) ? nothing : dB[i])
         for i in eachindex(specs)])
    src = reduce(vcat, [fill(speclabel(specs[i]), length(t[i])) for i in eachindex(specs)])
    ds = Dataset1D(Planes(traces, vars), defaultnoisecentre(first(specs)),
        speclabel(first(specs)), src)

    # named `regs`, the accessor `regions` being shadowed here by the keyword of that name
    regs = @something(regions, [defaultregion(ds)])
    expt = NutationExperiment(ds; phase, regions=regs)
    # `durations` applies to every spectrum, so it is only worth recording when the lists
    # agree; where they differ, each spectrum's own annotation reproduces the analysis.
    call = analysiscall("calibration1d", given;
        durations=(allequal(t) ? first(t) : nothing), phase, power=dB)
    return expt, ds, regs, call
end

"""
    calibrationanalysis(specs; kwargs...) -> (calibration, save)

Fit nutation calibration experiments without a window, returning the
[`B1Calibration`](@ref) they measure and a function that writes the analysis into a folder,
given one.

Whether and where those results are kept is for the analysis that asked for the calibration
to decide: a calibration fitted on the way to an exchange fit belongs inside that fit's
output folder, written when it is written, and not written at all if the fit is abandoned.
Keywords are [`calibration1d`](@ref)'s.

# Example
```julia
cal, save = calibrationanalysis(["cal/1", "cal/2"])
save(joinpath("out", "calibration"))     # fit.pdf, calibration.pdf, summary.txt, …
```
"""
function calibrationanalysis(specs::AbstractVector; kwargs...)
    expt, ds, regs, call = nutationexperiment(specs; kwargs...)
    results = analyse(expt)
    return (B1Calibration(results; source=join(string.(specs), ", ")),
        folder -> saveanalysis(expt, ds, results, regs, folder; call))
end

calibrationanalysis(spec; kwargs...) = calibrationanalysis([spec]; kwargs...)

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
`1/(4ν)`, and `σ` is the width of the B₁ distribution as a *percentage* of ν, which is what
damps the nutation: averaging `sin(2πνt)` over a Gaussian distribution of ν with mean ν̄ and
standard deviation σν̄ gives `sin(2πν̄t)·exp(−½(2πσν̄t)²)`, a Gaussian decay rather than an
exponential one. Fitting the envelope this way makes σ the same quantity
[`B1Distribution`](@ref) later samples, instead of an exponential rate needing conversion.

σ is fitted in percent rather than as a fraction so that every inhomogeneity this analysis
reports is in the unit it is quoted in, and incidentally so that the three parameters are of
comparable magnitude, which is what a Levenberg-Marquardt step prefers.

Relaxation during the pulse damps the nutation too, and is not separated here: σ is
therefore an upper bound, worst at the lowest power where the pulses are longest, which is
why a multi-power calibration reports the smallest of its estimates.

Only σ² enters, so the fitted sign is arbitrary and `postfit!` reports |σ|.
"""
function DampedSinusoidModel(; phase::Symbol=:sine)
    trig = phase === :cosine ? cos : sin
    est(x, y) = [maximum(abs.(y)), estimatefrequency(x, y, trig), 5.0]
    return CurveFitModel((x, p) -> @.(p[1] * trig(2π * p[2] * x) *
        exp(-0.5 * (2π * (p[3] / 100) * p[2] * x)^2)),
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

# Stored in the units `PARAM_UNITS` names for them: a 90° pulse reads naturally in µs,
# never in seconds, and an inhomogeneity in percent, which is what the model fits.
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
        # only σ² enters the model, so the fitted sign means nothing
        setpost!(r, seriesname(:sigma, s.group), abs(param(r, seriesname(:sigma, s.group))))
    end

    # The smallest estimate across power levels: every one is inflated by relaxation during
    # the pulse, least so at high power where the pulses are shortest.
    σ = [param(r, seriesname(:sigma, s.group)) for s in r.series]
    setpost!(r, :inhomogeneity, σ[argmin(Measurements.value.(σ))])
    # With one power level there is nothing to choose between, so `sigma` and
    # `inhomogeneity` are one number under two names: keep the one it is named for.
    length(r.series) == 1 && delete!(r.parameters, :sigma)

    # The calibration curve, where several power levels were measured. `nu1ref` is the
    # field at `powerref` and `linearity` the fitted exponent relative to the ideal power
    # law; with one power level they would only repeat `nu` and `power`.
    if hasvar(p, :power) && length(r.series) > 1
        cal = B1Calibration(r)
        setpost!(r, :powerref, db(refpower(cal)))
        setpost!(r, :nu1ref, ν1ref(cal))
        setpost!(r, :linearity, linearity(cal))
    end

    # One power level's results together, then the next, then the calibration:
    # `liftparameters!` inserted the fitted coefficients series by series and `setpost!`
    # appends, which would otherwise list every 90° pulse below every frequency.
    ordered = OrderedDict{Symbol,Any}()
    for s in r.series, key in (:nu, :pulse90, :sigma, :A)
        name = seriesname(key, s.group)
        haskey(r.parameters, name) && (ordered[name] = r.parameters[name])
    end
    for (name, value) in r.parameters
        haskey(ordered, name) || (ordered[name] = value)
    end
    r.parameters = ordered
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
`calibration1d`). To keep the analysis behind it as well as the calibration itself, use
[`calibrationanalysis`](@ref), which hands back both. Where several regions were
integrated, the first is the calibration; the others are presumably there for comparison.

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
    haskey(r.parameters, seriesname(:nu, first(r.series).group)) &&
        haskey(r.parameters, :inhomogeneity) ||
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
    # The headline inhomogeneity is already the smallest of the per-power estimates, and
    # `postfit!` reports it in the percent it is quoted in; a calibration holds a fraction.
    σ = Measurements.value(param(r, :inhomogeneity))
    return B1Calibration(power, ν1; nuc, inhomogeneity=σ / 100, source)
end

function B1Calibration(results::AbstractVector{RegionResult}; kwargs...)
    isempty(results) && throw(ArgumentError("no regions were analysed"))
    return B1Calibration(first(results); kwargs...)
end

function B1Calibration(specs::AbstractVector; kwargs...)
    return first(calibrationanalysis(specs; kwargs...))
end

function B1Calibration(spec::Union{AbstractString,Integer}; kwargs...)
    return B1Calibration([spec]; kwargs...)
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::NutationExperiment) = "Nutation calibration"

resultxfactor(e::NutationExperiment) = timescale(column(dataset(e).planes, :duration))[1]

function resultlabels(e::NutationExperiment)
    _, unit = timescale(column(dataset(e).planes, :duration))
    return ("Pulse duration / $unit", "Intensity")
end

function spectruminfo(::NutationExperiment, vars::NamedTuple)
    info = "$(round(1e6 * vars.duration; digits=1)) µs pulse"
    haskey(vars, :power) && return "$info at $(vars.power) dB"
    return info
end

"""
    resultfigure(expt::NutationExperiment, result, labels) -> Figure

One panel per power level rather than one axis for all of them.

The whole point of calibrating at several powers is that they nutate at different rates,
so their pulse durations differ by the same factor: 123 µs against 1220 µs here, which on
a shared axis leaves the fast series as a spike against the origin. Each series gets its
own axis, autoscaled to its own durations, stacked in the order they were measured and
labelled once underneath.
"""
function resultfigure(e::NutationExperiment, result, labels)
    xl, yl = resultlabels(e)
    fig = Figure()
    axes = Axis[]
    for label in labels, (i, s) in enumerate(resultplotdata(e, result, label))
        title = isempty(s.label) ? label : "$label ($(s.label))"
        ax = Axis(fig[length(axes)+1, 1]; ylabel=yl, title=title, titlesize=12,
            xgridvisible=false, ygridvisible=false)
        hlines!(ax, [0]; color=:grey)
        plotseries!(ax, s, seriescolor(i))
        push!(axes, ax)
    end
    # the shared quantity is labelled once, under the bottom panel
    isempty(axes) || (last(axes).xlabel = xl)
    return fig
end

"""
    calibrationplot(cal) -> Figure

The calibration curve: the field measured at each power level, the fitted power law through
them, and the residuals from it.

ν₁ is drawn on a log axis against power in dB, the coordinates in which the ideal law is a
straight line, so amplifier compression shows as curvature. The deviations are a percent or
two and invisible against a decade of field strength, hence the residual panel underneath,
in percent rather than in σ: what matters here is how far the amplifier departs from the
law, not how that compares with the fitting uncertainty.
"""
function calibrationplot(cal::B1Calibration)
    x = db.(powers(cal))
    y = Measurements.value.(fields(cal))
    yerr = Measurements.uncertainty.(fields(cal))
    fitted = [hz(p, cal) for p in powers(cal)]
    grid = range(minimum(x) - 0.5, maximum(x) + 0.5; length=100)

    fig = Figure()
    common = (; xgridvisible=false, ygridvisible=false)
    # data above, residuals below at a third the height, sharing an x-axis (matching the
    # Exchange1D result plots); the shared axis is labelled once, underneath
    ax1 = Axis(fig[1, 1]; ylabel="ν₁ / Hz", yscale=log10,
        title="Calibration: linearity $(fmt(linearity(cal))), " *
            "B₁ inhomogeneity $(round(100 * inhomogeneity(cal); digits=1))%",
        common...)
    ax2 = Axis(fig[2, 1]; xlabel="Power / dB", ylabel="Residual / %", common...)
    linkxaxes!(ax1, ax2)
    rowsize!(fig.layout, 1, Auto(false, 3))
    rowsize!(fig.layout, 2, Auto(false, 1))
    rowgap!(fig.layout, 4)

    lines!(ax1, grid, [hz(Power(g, :dB), cal) for g in grid]; color=:grey)
    errorbars!(ax1, x, y, yerr; whiskerwidth=8, color=seriescolor(1))
    scatter!(ax1, x, y; color=seriescolor(1))

    hlines!(ax2, [0]; color=:grey)
    errorbars!(ax2, x, 100 .* (y ./ fitted .- 1), 100 .* yerr ./ fitted;
        whiskerwidth=8, color=seriescolor(1))
    scatter!(ax2, x, 100 .* (y ./ fitted .- 1); color=seriescolor(1))
    return fig
end

"""
    saveextras!(expt::NutationExperiment, results, folder)

Write `calibration.pdf`, the calibration curve, where several power levels were measured.
One power level is a point rather than a curve, and gets no plot. As with
[`B1Calibration`](@ref), the first region carrying a curve is the calibration; any others
are there for comparison.
"""
function saveextras!(e::NutationExperiment, results, folder::AbstractString)
    hasvar(dataset(e).planes, :power) || return nothing
    i = findfirst(r -> length(r.series) > 1, results)
    isnothing(i) && return nothing
    path = joinpath(folder, "calibration.pdf")
    save(path, calibrationplot(B1Calibration(results[i])); backend=CairoMakie)
    return path
end

# Own display names and units, not the shared PARAM_LABELS/PARAM_UNITS tables -
# everything about this experiment's presentation lives here. The keys are ASCII (:nu, not
# :ν) because they become CSV column headers; the typeset names are in the labels.
#
# :sigma is the width of the B₁ distribution fitted at one power level and :inhomogeneity
# the smallest of those estimates, both in percent. They are named apart so that the
# headline figure is distinguishable from the estimate it was chosen from.
const NUTATION_PARAM_LABELS = Dict(:nu => "Nutation frequency",
    :pulse90 => "90° pulse",
    :sigma => "B₁ inhomogeneity",
    :inhomogeneity => "B₁ inhomogeneity",
    :power => "Power",
    :powerref => "Reference power",
    :nu1ref => "ν₁ at reference power",
    :linearity => "Linearity")
const NUTATION_PARAM_UNITS = Dict(:nu => "Hz",
    :pulse90 => "us",
    :sigma => "%",
    :inhomogeneity => "%",
    :power => "dB",
    :powerref => "dB",
    :nu1ref => "Hz")

function paramlabel(e::NutationExperiment, name::Symbol)
    # Several power levels each report a B₁ inhomogeneity, and the region-level figure is
    # the smallest of them: worth saying, where a single power level has nothing to choose
    # between and the plain name is right.
    name === :inhomogeneity && !isempty(groupcols(e)) && return "B₁ inhom. (smallest)"
    return get(NUTATION_PARAM_LABELS, name, get(PARAM_LABELS, name, string(name)))
end
function paramunit(::NutationExperiment, name::Symbol)
    return get(NUTATION_PARAM_UNITS, name, get(PARAM_UNITS, name, ""))
end
