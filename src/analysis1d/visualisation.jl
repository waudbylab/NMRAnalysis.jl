# Result-panel visualisation. Pure data builders (result -> plot primitives) dispatched on
# the experiment type, plus axis labels and a summary string. The GUI lifts Observables off
# these; saving reuses the same builders, so live and exported plots share one code path.

"""
    ResultSeries

One coloured curve in the result panel: observed `points` (`Vector{Point2f}`), `errors`
(`Vector{NTuple{3,Float64}}` of `(x, y, σ)`), a `fitline` (`Vector{Point2f}`, empty if
unfitted), and a `label` (e.g. the TRACT `which`, or the STD saturation frequency; empty
for an ungrouped single series).
"""
struct ResultSeries
    points::Vector{Point2f}
    errors::Vector{Tuple{Float64,Float64,Float64}}
    fitline::Vector{Point2f}
    label::String
end

"""
    timescale(times) -> (factor, unit)

Pick a display unit (s / ms / µs) and its multiplying factor for a set of raw times (in
s), so values don't show as awkward tiny fractions - e.g. nutation's calibration pulses
(tens of µs to ~1 ms) read far more naturally in µs than in s, just as sub-2s relaxation
delays read better in ms.
"""
function timescale(times)
    isempty(times) && return (1.0, "s")
    m = maximum(times)
    m < 2e-3 && return (1.0e6, "µs")
    m < 2.0 && return (1.0e3, "ms")
    return (1.0, "s")
end

"""
    resultxfactor(expt) -> Float64

Multiplier applied to the result panel's x-values (and its axis label, via
`result_labels`) for display - `1.0` (unscaled) for every experiment except those with a
time-valued x-axis (relaxation delays, TRACT delays, nutation pulse durations), which
switch to whichever of s/ms/µs suits their actual scale (see `timescale`). Curve-fitting
itself always uses the raw, unscaled values; this only affects what's plotted/exported.
"""
resultxfactor(::Experiment1D) = 1.0
resultxfactor(e::RelaxationExperiment) = timescale(column(dataset(e).planes, :time))[1]
resultxfactor(e::TractExperiment) = timescale(column(dataset(e).planes, :time))[1]
resultxfactor(e::NutationExperiment) = timescale(column(dataset(e).planes, :duration))[1]

"""
    result_plotdata(expt, result, activelabel) -> Vector{ResultSeries}

Plot primitives for the active region, one `ResultSeries` per group (e.g. TRACT's
TROSY/anti-TROSY pair, or STD's saturation frequencies) so each can be drawn in a
distinct, matching colour. The generic method covers the curve-fit / NoFitting
experiments; STD overrides it.
"""
function result_plotdata(e::Experiment1D, result, activelabel::AbstractString)
    series = filter(s -> s.region == activelabel, result.series)
    factor = resultxfactor(e)
    return map(series) do s
        points = Point2f.(factor .* s.x, Measurements.value.(s.y))
        errors = [(factor * s.x[k], Measurements.value(s.y[k]), Measurements.uncertainty(s.y[k]))
                  for k in eachindex(s.x)]
        fitline = if !(s.model isa NoFitting) && !isempty(s.params)
            xs = collect(range(min(0.0, minimum(s.x)), 1.05 * maximum(s.x), 100))
            Point2f.(factor .* xs, s.model.func(xs, Measurements.value.(s.params)))
        else
            Point2f[]
        end
        ResultSeries(points, errors, fitline, groupname(s.group))
    end
end

function result_plotdata(::STDExperiment, result, activelabel::AbstractString)
    pts = filter(p -> p.region == activelabel, result.points)
    sats = unique(p.sat for p in pts)
    return map(sats) do sat
        ps = filter(p -> p.sat == sat, pts)
        points = [Point2f(p.tsat, Measurements.value(p.std)) for p in ps]
        errors = [(p.tsat, Measurements.value(p.std), Measurements.uncertainty(p.std))
                  for p in ps]
        bi = findfirst(b -> b.region == activelabel && b.sat == sat, result.buildups)
        fitline = if !isnothing(bi)
            b = result.buildups[bi]
            tmax = maximum(p.tsat for p in ps)
            xs = collect(range(0.0, 1.05 * tmax, 100))
            smax, k = Measurements.value(b.std_af_max), Measurements.value(b.k)
            Point2f.(xs, @.(smax * (1 - exp(-k * xs))))
        else
            Point2f[]
        end
        ResultSeries(points, errors, fitline, string(sat))
    end
end

"""Human-readable label for a grouping key, e.g. `(which = :trosy,)` → `"trosy"`."""
function groupname(group::NamedTuple)
    isempty(group) && return ""
    return join((string(v) for v in values(group)), ", ")
end

# axis labels for the result panel
result_labels(::Experiment1D) = ("x", "Integrated intensity (a.u.)")
function result_labels(e::RelaxationExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end
function result_labels(e::TractExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end
function result_labels(e::NutationExperiment)
    _, unit = timescale(column(dataset(e).planes, :duration))
    return ("Pulse duration / $unit", "Integrated intensity (a.u.)")
end
result_labels(::KineticsExperiment) = ("Time", "Integrated intensity (a.u.)")
result_labels(::DiffusionExperiment) = ("Relative gradient strength", "Integrated intensity (a.u.)")
result_labels(::STDExperiment) = ("Saturation time / s", "STD fraction")

"""
    seriesnames(expt) -> Union{Nothing,Vector{String}}

Fixed group names (in `seriescolor` order) for the result panel's legend, when the
experiment's groups are few and consistently named (e.g. TRACT's TROSY/anti-TROSY
pair) - `axislegend` is used instead of labelling each curve directly. `nothing` (the
default) keeps the inline curve labels, which suit experiments whose group count/names
vary per dataset (e.g. STD's saturation frequencies).
"""
seriesnames(::Experiment1D) = nothing
seriesnames(::TractExperiment) = ["TROSY", "anti-TROSY"]

# analysis-type name shown in the GUI's bold in-panel title, e.g. "NMRAnalysis: Relaxation"
windowtitle(::Experiment1D) = "1D analysis"
windowtitle(::RelaxationExperiment) = "Relaxation"
windowtitle(::TractExperiment) = "TRACT"
windowtitle(::NutationExperiment) = "Nutation calibration"
windowtitle(::KineticsExperiment) = "Kinetics"
windowtitle(::DiffusionExperiment) = "Diffusion"
windowtitle(::STDExperiment) = "STD"

"""
    spectruminfo(expt, vars) -> String

One-line description of a single plane, built from its arrayed variables `vars` (the
`NamedTuple` for that plane, e.g. `(; time=0.4, which=:anti)`) - shown next to the
spectrum slider so the displayed slice is identifiable at a glance (e.g. "0.4 s delay
(anti-TROSY)"). The generic fallback lists every variable as `key=value`; each
experiment overrides it with wording natural to its own variables.
"""
function spectruminfo(::Experiment1D, vars::NamedTuple)
    isempty(vars) && return ""
    return join(("$k=$v" for (k, v) in pairs(vars)), ", ")
end

spectruminfo(::RelaxationExperiment, vars::NamedTuple) = "$(round(vars.time; digits=3)) s delay"

function spectruminfo(::TractExperiment, vars::NamedTuple)
    which = vars.which == :trosy ? "TROSY" : "anti-TROSY"
    return "$(round(vars.time; digits=3)) s delay ($which)"
end

function spectruminfo(::NutationExperiment, vars::NamedTuple)
    return "$(round(1e6 * vars.duration; digits=1)) µs pulse"
end

function spectruminfo(::KineticsExperiment, vars::NamedTuple)
    s = "$(round(vars.time; digits=3)) s"
    haskey(vars, :run) && (s *= " (run $(vars.run))")
    return s
end

function spectruminfo(::DiffusionExperiment, vars::NamedTuple)
    return "gradient = $(round(100 * vars.gradient; digits=1))%"
end

spectruminfo(::STDExperiment, vars::NamedTuple) = "$(vars.sat), $(round(vars.tsat; digits=3)) s sat"

# ---- summary text -------------------------------------------------------------

"""Compact, consistently-rounded rendering of a `Measurement` (or plain number)."""
fmt(x::Measurement, digits=4) = string(round(Measurements.value(x); sigdigits=digits),
                                       " ± ",
                                       round(Measurements.uncertainty(x); sigdigits=2))
fmt(x::Real, digits=4) = string(round(x; sigdigits=digits))
fmt(::Nothing, digits=4) = "n/a"

const PARAM_UNITS = Dict("R" => " s⁻¹", "ν" => " Hz")

"""
    resultsheader(expt, result, activelabel) -> String

Just the active region's raw per-series fitted parameters (e.g. `A`, `R`) - the
"primary" half of `summary_text`, split out so the GUI can show it and the "second-stage"
derived summary (see `secondarytext`, e.g. TRACT's τc) as visually separate blocks.
Doesn't repeat the region name, group headers only (e.g. TRACT's "trosy"/"anti"), since
the region itself is already named elsewhere in the GUI (the fit axis title).
"""
function resultsheader(e::Experiment1D, result, activelabel::AbstractString)
    io = IOBuffer()
    for s in result.series
        s.region == activelabel || continue
        isempty(s.group) || println(io, groupname(s.group))
        for (name, p) in zip(s.names, s.params)
            println(io, "  $(rpad(name, 6)) = $(fmt(p))$(get(PARAM_UNITS, name, ""))")
        end
        println(io)
    end
    return String(take!(io))
end
resultsheader(::STDExperiment, result, activelabel::AbstractString) = ""

"""
    summary_text(expt, result, activelabel) -> String

Human-readable summary of the fit for the active region only, shown in the GUI and
written to `summary.txt` (matching the result panel and plot, which are likewise
restricted to the active region).

# Note
Output-format consistency across the different 1D/2D analyses (units, significant
figures, CSV vs text) is an open question — see `PLAN.md` — this only tidies the
single-experiment summary shown here.
"""
function summary_text(e::Experiment1D, result, activelabel::AbstractString)
    io = IOBuffer()
    for s in result.series
        s.region == activelabel || continue
        header = isempty(s.group) ? s.region : "$(s.region) ($(groupname(s.group)))"
        println(io, header)
        println(io, "-"^length(header))
        for (name, p) in zip(s.names, s.params)
            println(io, "  $(rpad(name, 6)) = $(fmt(p))$(get(PARAM_UNITS, name, ""))")
        end
        println(io)
    end
    _summary_extra(io, e, result.summary, activelabel)
    return String(take!(io))
end

_summary_extra(::IOBuffer, ::Experiment1D, ::Nothing, ::AbstractString) = nothing
# Same "  name  = value" shape as `resultsheader`'s per-series parameters, and no region
# name repeated - the region is already named once, above both blocks, in the GUI (and
# in `summary_text`'s own per-series header, when saved to a file).
function _summary_extra(io::IOBuffer, ::TractExperiment, summary, activelabel)
    for s in summary
        s.region == activelabel || continue
        println(io, "  $(rpad("η", 6)) = $(fmt(s.ηxy)) s⁻¹")
        println(io, "  $(rpad("τc", 6)) = $(fmt(s.τc)) ns")
    end
end
function _summary_extra(io::IOBuffer, ::NutationExperiment, summary, activelabel)
    for s in summary
        s.region == activelabel || continue
        println(io, "  $(rpad("ν₁", 6)) = $(fmt(s.ν)) Hz")
        println(io, "  $(rpad("90°", 6)) = $(fmt(1e6 * s.pulse90)) µs")
    end
end
function _summary_extra(io::IOBuffer, ::DiffusionExperiment, summary, activelabel)
    for s in summary
        s.region == activelabel || continue
        println(io, "  $(rpad("D", 6)) = $(fmt(s.D)) m² s⁻¹")
        isnothing(s.rH) || println(io, "  $(rpad("rH", 6)) = $(fmt(s.rH)) Å")
        isnothing(s.η) || println(io, "  $(rpad("η", 6)) = $(fmt(s.η)) mPa s")
    end
end

"""
    secondarytext(expt, result, activelabel) -> String

Formatted text for the active region's "second-stage" derived summary alone (e.g.
TRACT's τc and ΔR, computed by combining its trosy/anti series) - empty when the
experiment has none (`result.summary === nothing`, e.g. relaxation/kinetics) or the
active region isn't ready yet (e.g. only one of a TRACT pair fitted so far). Reuses
`_summary_extra`, which already knows how to format and filter this per experiment type.
"""
function secondarytext(e::Experiment1D, result, activelabel::AbstractString)
    io = IOBuffer()
    _summary_extra(io, e, result.summary, activelabel)
    return String(take!(io))
end
# STD's result has no `.series`/`.summary` fields to split into a resultsheader/
# secondarytext pair (it's a contrast rather than a curve-fit experiment - see
# `resultsheader`'s STD override above, near `PARAM_UNITS`) - show its whole summary
# here instead.
secondarytext(e::STDExperiment, result, activelabel::AbstractString) = summary_text(e, result,
                                                                                    activelabel)

function summary_text(::STDExperiment, result, activelabel::AbstractString)
    io = IOBuffer()
    println(io, "STD fractions")
    println(io, "-------------")
    for p in result.points
        p.region == activelabel || continue
        println(io, "  $(p.region) / $(p.sat) @ $(p.tsat) s : $(fmt(p.std))")
    end
    buildups = filter(b -> b.region == activelabel, result.buildups)
    if !isempty(buildups)
        println(io, "\nBuildup (initial slope, T1-corrected)")
        println(io, "--------------------------------------")
        for b in buildups
            println(io,
                    "  $(b.region) / $(b.sat) : STD-AF₀ = $(fmt(b.std_af0)), k = $(fmt(b.k)) s⁻¹")
        end
    end
    epitope = filter(ep -> ep.region == activelabel, result.epitope)
    if !isempty(epitope)
        println(io, "\nEpitope map (relative to strongest signal)")
        println(io, "-------------------------------------------")
        for ep in epitope
            println(io, "  $(ep.region) / $(ep.sat) : $(round(100 * ep.relative; digits=0)) %")
        end
    end
    return String(take!(io))
end
