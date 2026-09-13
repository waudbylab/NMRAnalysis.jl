# Result-panel visualisation: the generic half. Pure data builders (result -> plot
# primitives), axis labels, and the summary formatter, each with a default that covers
# the curve-fit experiments. Per-experiment overrides live in the `expt-*.jl` files.
#
# The GUI lifts Observables off these builders; saving reuses them, so live and exported
# plots share one code path.

"""
    ResultSeries

One coloured curve in the result panel: observed `points` (`Vector{Point2f}`), `errors`
(`Vector{NTuple{3,Float64}}` of `(x, y, σ)`), a `fitline` (`Vector{Point2f}`, empty if
unfitted), and a `label` (e.g. the TRACT `which`, or kinetics' `run`; empty for an
ungrouped single series).
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

Multiplier applied to the result panel's x-values for display, `1.0` except where a
time-valued axis reads better in ms or µs (see [`timescale`](@ref)). Fitting always uses
the raw values.
"""
resultxfactor(::Experiment1D) = 1.0

"""
    resultplotdata(expt, result, activelabel) -> Vector{ResultSeries}

Plot primitives for the active region, one `ResultSeries` per series, so each is drawn in
its own colour.
"""
function resultplotdata(e::Experiment1D, result, activelabel::AbstractString)
    i = findfirst(r -> r.region == activelabel, result)
    isnothing(i) && return ResultSeries[]
    factor = resultxfactor(e)
    return map(result[i].series) do s
        points = Point2f.(factor .* s.x, Measurements.value.(s.y))
        errors = [(factor * s.x[k], Measurements.value(s.y[k]),
                   Measurements.uncertainty(s.y[k]))
                  for k in eachindex(s.x)]
        fitline = if !(s.model isa NoFitting) && !isempty(s.coefficients)
            xs = collect(range(min(0.0, minimum(s.x)), 1.05 * maximum(s.x), 100))
            Point2f.(factor .* xs, s.model.func(xs, Measurements.value.(s.coefficients)))
        else
            Point2f[]
        end
        return ResultSeries(points, errors, fitline, groupname(s.group))
    end
end

"""Human-readable label for a grouping key, e.g. `(which = :trosy,)` → `"trosy"`."""
function groupname(group::NamedTuple)
    isempty(group) && return ""
    return join((string(v) for v in values(group)), ", ")
end

"""
    show(io, result::RegionResult)

One line per region: its label followed by every parameter reported for it, so that
closing an analysis window prints a short table rather than pages of measured points.

Units are absent, being the one thing a `RegionResult` cannot know: `paramunit` is
dispatched on the experiment. The results panel, `summary.txt` and the CSVs carry them.
"""
function Base.show(io::IO, r::RegionResult)
    print(io, r.region)
    print(io, isconverged(r) ? ": " : " (not converged): ")
    isempty(r.parameters) &&
        return print(io, "$(sum(length(s.x) for s in r.series; init=0)) points, unfitted")
    return print(io, join(("$name = $value" for (name, value) in r.parameters), ", "))
end

"axis labels for the result panel"
resultlabels(::Experiment1D) = ("x", "Integrated intensity (a.u.)")

"""
    seriesnames(expt) -> Union{Nothing,Vector{String}}

Fixed group names (in `seriescolor` order) for the result panel's legend, when the
experiment's groups are few and consistently named (e.g. TRACT's TROSY/anti-TROSY
pair) - `axislegend` is used instead of labelling each curve directly. `nothing` (the
default) keeps the inline curve labels, which suit experiments whose group count/names
vary per dataset (e.g. kinetics' runs, when several are present).
"""
seriesnames(::Experiment1D) = nothing

"analysis-type name shown in the GUI's bold in-panel title, e.g. \"NMRAnalysis: Relaxation\""
windowtitle(::Experiment1D) = "1D analysis"

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

# ---- visualisation strategy ---------------------------------------------------

"""
    ResultVisualisation

How an experiment's results are drawn in the fit panel. A hierarchy orthogonal to
`Experiment1D`, joined by [`visualisationtype`](@ref), as GUI2D's `VisualisationStrategy`
is to its `Experiment`, so a presentation can be shared or swapped without touching the
science.

A strategy implements three methods:

- `completeresultstate!(state, expt, ::V)` — build the Observables the panel reads
- `resultpanel!(gui, state, expt, ::V)`    — the live GUI panel
- `plotresult!(ax, expt, result, label, i0, ::V)` — the static CairoMakie export

The live and export paths share one data getter (here `resultplotdata`) so that what is
saved is what was on screen.

[`SeriesVisualisation`](@ref) covers every curve-fit experiment and is the default.
"""
abstract type ResultVisualisation end

"""
    SeriesVisualisation()

Observed points with error bars plus a fitted line, one colour per group. The default,
and what relaxation, TRACT, nutation, diffusion and kinetics all use.
"""
struct SeriesVisualisation <: ResultVisualisation end

"""
    visualisationtype(expt) -> ResultVisualisation

The strategy drawing `expt`'s results. Defaults to [`SeriesVisualisation`](@ref), which
suits every curve-fit experiment; override it in an `expt-*.jl` only where that shape
genuinely will not do.
"""
visualisationtype(::Experiment1D) = SeriesVisualisation()

# Trait forwarding: callers use the three-argument forms and never name a strategy.
function completeresultstate!(state, expt::Experiment1D)
    return completeresultstate!(state, expt, visualisationtype(expt))
end
function resultpanel!(gui, state, expt::Experiment1D)
    return resultpanel!(gui, state, expt,
                        visualisationtype(expt))
end
function plotresult!(ax, expt::Experiment1D, result, label, i0=0)
    return plotresult!(ax, expt, result, label, i0, visualisationtype(expt))
end

# ---- colours ------------------------------------------------------------------

# Explicit RGBAf conversion, so every colour flowing into the same Observable/plot-array
# is the same concrete type as `Makie.wong_colors()` returns colours without alpha.
const PALETTE = [RGBAf(c.r, c.g, c.b, 1.0) for c in Makie.wong_colors()]

"""Colour for the `i`-th result series, cycling through a fixed palette."""
seriescolor(i) = PALETTE[mod1(i, length(PALETTE))]

# ---- summary text -------------------------------------------------------------

"""Compact, consistently-rounded rendering of a `Measurement` (or plain number)."""
function fmt(x::Measurement, digits=4)
    return string(round(Measurements.value(x); sigdigits=digits),
                  " ± ",
                  round(Measurements.uncertainty(x); sigdigits=2))
end
fmt(x::Real, digits=4) = string(round(x; sigdigits=digits))
fmt(::Nothing, digits=4) = "n/a"

# Only the parameters genuinely shared by more than one experiment. Anything specific to
# one lives in its own expt-*.jl, beside the `postfit!` that produces it.
#
# `:R` has a unit here but no label: every model producing a rate means s⁻¹ by it, but
# relaxation and TRACT mean a relaxation rate where nutation's damped sinusoid means a
# decay rate, so "Relaxation rate" would mislabel the latter.
const PARAM_UNITS = Dict(:R => "s-1")
const PARAM_LABELS = Dict(:A => "Amplitude")

"""
    paramlabel(expt, name) -> String

Display name for parameter `name` (fitted or derived). Experiment-dispatched: the default
falls back to the small shared table above, then to the bare symbol. An experiment that
introduces its own parameters overrides this with its own local table - see e.g.
`NUTATION_PARAM_LABELS` in `expt-nutation.jl`.
"""
paramlabel(::Experiment1D, name::Symbol) = get(PARAM_LABELS, name, string(name))

"""
    paramunit(expt, name) -> String

Unit for parameter `name`, dispatched and overridden as [`paramlabel`](@ref) is.

Units are held in **ASCII** (`"s-1"`, `"us"`, `"1e-10 m2/s"`), because that is the form the
CSV column headers carry and a file header must survive being opened on any machine - see
`docs/src/advanced/conventions.md`. [`prettyunit`](@ref) turns one into the typeset form
for the GUI panel and `summary.txt`, so there is one table rather than two that can drift.
"""
paramunit(::Experiment1D, name::Symbol) = get(PARAM_UNITS, name, "")

"""
    displaylabel(expt, name) -> String

[`paramlabel`](@ref) for a parameter as stored on a region, which for a fitted one carries
the series it came from: `:R_trosy` shows as "Relaxation rate (trosy)". The quantity comes
from [`baseparam`](@ref), so an experiment's label table needs only bare names.
"""
function displaylabel(e::Experiment1D, name::Symbol)
    base = baseparam(name)
    base === name && return paramlabel(e, name)
    return "$(paramlabel(e, base)) ($(string(name)[(length(string(base)) + 2):end]))"
end

"""
    prettyunit(unit) -> String

The typeset form of an ASCII unit, for display only: `"s-1"` → `" s⁻¹"`, `"us"` → `" µs"`.
Includes the leading space so it appends to a formatted value, and returns `""` for a
dimensionless quantity. Anything not in the table is shown as written (`"Hz"`, `"ns"`).
"""
function prettyunit(unit::AbstractString)
    isempty(unit) && return ""
    return " " * get(PRETTY_UNITS, unit, unit)
end

const PRETTY_UNITS = Dict("s-1" => "s⁻¹",
                          "us" => "µs",
                          "A" => "Å",
                          "m2/s" => "m² s⁻¹",
                          "1e-10 m2/s" => "×10⁻¹⁰ m² s⁻¹")

"Display unit for a parameter stored on a region: its quantity's unit, typeset."
prettyparamunit(e::Experiment1D, name::Symbol) = prettyunit(paramunit(e, baseparam(name)))

"""
    coordinateunit(expt, name) -> String

ASCII unit for a coordinate, the fit axis or a grouping variable, as it appears in a column
header. Kept apart from [`paramunit`](@ref) because the two name spaces are independent:
nothing stops an experiment arraying a variable named like a fitted parameter. Categorical
coordinates (`:which`, `:run`) and relative ones (`:gradient`) have no unit.
"""
coordinateunit(::Experiment1D, name::Symbol) = get(COORDINATE_UNITS, name, "")

const COORDINATE_UNITS = Dict(:time => "s", :duration => "s")

"""
    paramblock(io, expt, params, width=nothing)

Write every entry of a parameter dictionary as `name value unit`, in insertion order, the
names padded to `width` (computed from `params` alone when omitted) so the values line up.
The GUI panel passes a `width` computed across every region ([`panelwidth`](@ref)) so the
column does not jump as regions are selected; `summary.txt` leaves it to align to its own
block.
"""
function paramblock(io::IO, expt::Experiment1D, params, width=nothing)
    isempty(params) && return nothing
    w = something(width,
                  maximum(length(displaylabel(expt, name)) for name in keys(params)) + 2)
    for (name, value) in params
        println(io,
                "$(rpad(displaylabel(expt, name), w))$(fmt(value))$(prettyparamunit(expt, name))")
    end
    return nothing
end

"""
    paramtext(expt, params, width=nothing) -> String

`paramblock`, returned as a string rather than written to an `IO` - the form the GUI
panel wants, built from the same alignment logic `summarytext` uses for `summary.txt` so
the two never drift apart. Empty (not a blank line) when `params` is empty, so a caller
can skip a block that has nothing to show.
"""
function paramtext(expt::Experiment1D, params, width=nothing)
    io = IOBuffer()
    paramblock(io, expt, params, width)
    return String(take!(io))
end

"""
    panelwidth(expt, result, activelabel) -> Int

The label-column width for the results panel: the longest display label among the
parameters shown for `activelabel`, plus a gap, or `nothing` when there is nothing to show
yet. Computed across the whole panel so "Amplitude (trosy)" and "Correlation time (τc)"
share one column.
"""
function panelwidth(expt::Experiment1D, result, activelabel::AbstractString)
    len = 0
    for r in result
        r.region == activelabel || continue
        for k in keys(r.parameters)
            len = max(len, length(displaylabel(expt, k)))
        end
    end
    return len == 0 ? nothing : len + 2
end

"""
    plaintext(text) -> RichText

Wrap `text` in an explicit `font=:regular` span. A `RichText` font is not scoped to each
sibling: a plain `String` child inherits whatever font the previous sibling left active, so
without this a parameter block following a bold heading renders bold too.
"""
plaintext(text::AbstractString) = rich(text; font=:regular)

"""
    BLANK_RICHTEXT

Placeholder for "nothing to show". An empty `RichText` renders zero glyphs, and Makie's
`GlyphCollection` cannot build itself from a zero-length glyph vector: it cannot infer
`rotations` as `Vector{Quaternionf}` from an empty comprehension, and errors rather than
rendering blank. A single space has one glyph and is invisible. Without it, `resultstext`
crashes for every region of a `NoFitting` experiment and any region not yet fitted.
"""
const BLANK_RICHTEXT = rich(" ")

"""
    richtext(spans) -> RichText

`rich(spans...)`, substituting [`BLANK_RICHTEXT`](@ref) when `spans` is empty - see its
docstring for why that substitution is load-bearing, not decorative.
"""
richtext(spans) = isempty(spans) ? BLANK_RICHTEXT : rich(spans...)

"""
    resultstext(expt, result, activelabel, width=nothing) -> RichText

Every parameter reported for the active region, as one block - the fitted ones under the
names their series gave them, then whatever the experiment derived from them. Doesn't
repeat the region name, which is already shown elsewhere in the panel.

`width` aligns the values to a fixed column - see [`panelwidth`](@ref); omitted, the block
aligns only to its own labels.
"""
function resultstext(e::Experiment1D, result, activelabel::AbstractString, width=nothing)
    spans = Any[]
    for r in result
        r.region == activelabel || continue
        block = paramtext(e, r.parameters, width)
        isempty(block) && continue
        push!(spans, plaintext(block), "\n")
    end
    return richtext(spans)
end

"""
    summarytext(expt, result, activelabel) -> String

Human-readable summary of the fit for the active region only, shown in the GUI and
written to `summary.txt` (matching the result panel and plot, which are likewise
restricted to the active region).

# Note
Output-format consistency across the different 1D/2D analyses (units, significant
figures, CSV vs text) is an open question — see `PLAN.md`.
"""
function summarytext(e::Experiment1D, result, activelabel::AbstractString)
    io = IOBuffer()
    for r in result
        r.region == activelabel || continue
        println(io, r.region)
        println(io, "-"^length(r.region))
        paramblock(io, e, r.parameters)
        println(io)
    end
    return String(take!(io))
end

# ---- SeriesVisualisation ------------------------------------------------------

"""
    completeresultstate!(state, expt, ::SeriesVisualisation)

Per-series plot data for the active region (one entry per group, e.g. TRACT's
TROSY/anti-TROSY), flattened into the contiguous, coloured arrays Makie plot objects want -
one `scatter!`/`lines!` per panel rather than a variable number of separate plot calls.
"""
function completeresultstate!(state, expt::Experiment1D, ::SeriesVisualisation)
    state[:seriesdata] = lift(state[:result], state[:activelabel]) do res, lbl
        return resultplotdata(expt, res, lbl)
    end
    state[:flatpoints] = lift(sd -> reduce(vcat, (s.points for s in sd); init=Point2f[]),
                              state[:seriesdata])
    state[:flatpointcolors] = lift(state[:seriesdata]) do sd
        return reduce(vcat,
                      (fill(seriescolor(i), length(s.points)) for (i, s) in enumerate(sd));
                      init=RGBAf[])
    end
    state[:flaterrors] = lift(state[:seriesdata]) do sd
        return reduce(vcat, (s.errors for s in sd);
                      init=Tuple{Float64,Float64,Float64}[])
    end
    state[:flaterrorcolors] = state[:flatpointcolors]
    state[:flatfit] = lift(state[:seriesdata], state[:isfitting]) do sd, fitting
        fitting || return Point2f[]
        pts = Point2f[]
        for s in sd
            isempty(s.fitline) && continue
            append!(pts, s.fitline)
            push!(pts, Point2f(NaN32, NaN32))
        end
        return pts
    end
    state[:flatfitcolors] = lift(state[:seriesdata], state[:isfitting]) do sd, fitting
        fitting || return RGBAf[]
        cols = RGBAf[]
        for (i, s) in enumerate(sd)
            isempty(s.fitline) && continue
            append!(cols, fill(seriescolor(i), length(s.fitline) + 1))
        end
        return cols
    end
    state[:seriestextpos] = lift(sd -> [isempty(s.points) ? Point2f(NaN32, NaN32) :
                                        last(s.points) for s in sd], state[:seriesdata])
    state[:seriestexttxt] = lift(sd -> [s.label for s in sd], state[:seriesdata])
    state[:seriestextcolor] = lift(sd -> [seriescolor(i) for i in eachindex(sd)],
                                   state[:seriesdata])
    return state
end

"""
    resultpanel!(gui, state, expt, ::SeriesVisualisation)

Build the live fit panel into `gui[:panelresult]`, and record its axis as `gui[:axfit]`.
"""
function resultpanel!(gui, state, expt::Experiment1D, ::SeriesVisualisation)
    xl, yl = resultlabels(expt)
    fittitle = lift(lbl -> isempty(lbl) ? "Fit" : "Fit: $lbl", state[:activelabel])
    ax = Axis(gui[:panelresult]; xlabel=xl, ylabel=yl, title=fittitle,
              xgridvisible=false, ygridvisible=false)
    gui[:axfit] = ax
    hlines!(ax, [0]; color=:grey)
    errorbars!(ax, state[:flaterrors]; whiskerwidth=8, color=state[:flaterrorcolors])
    scatter!(ax, state[:flatpoints]; color=state[:flatpointcolors])
    lines!(ax, state[:flatfit]; color=state[:flatfitcolors])
    # Most experiments label each curve directly (group counts/names vary, e.g.
    # kinetics' runs); experiments with a fixed, small set of named groups (e.g. TRACT's
    # TROSY/anti-TROSY) get a proper axislegend instead via `seriesnames`.
    legendnames = seriesnames(expt)
    if isnothing(legendnames)
        text!(ax, state[:seriestextpos]; text=state[:seriestexttxt],
              color=state[:seriestextcolor], fontsize=11, align=(:left, :center),
              offset=(6, 0))
    else
        legendelements = [MarkerElement(; color=seriescolor(i), marker=:circle)
                          for i in eachindex(legendnames)]
        axislegend(ax, legendelements, legendnames; position=:rt)
    end
    on(_ -> autolimits!(ax), state[:seriesdata])
    return ax
end

"""
    plotresult!(ax, expt, result, label, i0, ::SeriesVisualisation) -> Int

Draw region `label` into a static axis, starting colour indices at `i0 + 1`; returns the
number of series drawn, so callers can decide whether a legend is worthwhile. Shares
`resultplotdata` with the live panel, so an exported plot matches what was on screen.
"""
function plotresult!(ax, expt::Experiment1D, result, label, i0, ::SeriesVisualisation)
    i = i0
    for s in resultplotdata(expt, result, label)
        i += 1
        c = seriescolor(i)
        lbl = isempty(s.label) ? label : "$label ($(s.label))"
        errorbars!(ax, s.errors; whiskerwidth=8, color=c)
        scatter!(ax, s.points; color=c, label=lbl)
        isempty(s.fitline) || lines!(ax, s.fitline; color=c)
    end
    return i - i0
end
