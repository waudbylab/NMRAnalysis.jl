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
`resultlabels`) for display - `1.0` (unscaled) for every experiment except those with a
time-valued x-axis (relaxation delays, TRACT delays, nutation pulse durations), which
switch to whichever of s/ms/µs suits their actual scale (see `timescale`). Curve-fitting
itself always uses the raw, unscaled values; this only affects what's plotted/exported.
"""
resultxfactor(::Experiment1D) = 1.0

"""
    resultplotdata(expt, result, activelabel) -> Vector{ResultSeries}

Plot primitives for the active region, one `ResultSeries` per group (e.g. TRACT's
TROSY/anti-TROSY pair, or STD's saturation frequencies) so each can be drawn in a
distinct, matching colour. The generic method covers the curve-fit / NoFitting
experiments; STD overrides it.
"""
function resultplotdata(e::Experiment1D, result, activelabel::AbstractString)
    series = filter(s -> s.region == activelabel, result)
    factor = resultxfactor(e)
    return map(series) do s
        points = Point2f.(factor .* s.x, Measurements.value.(s.y))
        errors = [(factor * s.x[k], Measurements.value(s.y[k]), Measurements.uncertainty(s.y[k]))
                  for k in eachindex(s.x)]
        fitline = if !(s.model isa NoFitting) && !isempty(s.parameters)
            xs = collect(range(min(0.0, minimum(s.x)), 1.05 * maximum(s.x), 100))
            ps = Measurements.value.(collect(values(s.parameters)))
            Point2f.(factor .* xs, s.model.func(xs, ps))
        else
            Point2f[]
        end
        ResultSeries(points, errors, fitline, groupname(s.group))
    end
end

"""Human-readable label for a grouping key, e.g. `(which = :trosy,)` → `"trosy"`."""
function groupname(group::NamedTuple)
    isempty(group) && return ""
    return join((string(v) for v in values(group)), ", ")
end

"axis labels for the result panel"
resultlabels(::Experiment1D) = ("x", "Integrated intensity (a.u.)")

"""
    seriesnames(expt) -> Union{Nothing,Vector{String}}

Fixed group names (in `seriescolor` order) for the result panel's legend, when the
experiment's groups are few and consistently named (e.g. TRACT's TROSY/anti-TROSY
pair) - `axislegend` is used instead of labelling each curve directly. `nothing` (the
default) keeps the inline curve labels, which suit experiments whose group count/names
vary per dataset (e.g. STD's saturation frequencies).
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

How an experiment's results are drawn in the fit panel. A hierarchy deliberately
*orthogonal* to `Experiment1D`, joined by [`visualisationtype`](@ref), exactly as GUI2D's
`VisualisationStrategy` is orthogonal to its `Experiment` - so that a presentation can be
shared by unrelated experiments, or swapped without touching the science.

A strategy implements three methods, the same trio GUI2D uses:

- `completeresultstate!(state, expt, ::V)` — build the Observables the panel reads
- `resultpanel!(gui, state, expt, ::V)`    — the live GUI panel
- `plotresult!(ax, expt, result, label, i0, ::V)` — the static CairoMakie export

The live and export paths share one data getter (here `resultplotdata`) so that what is
saved is what was on screen. A strategy's getter may return whatever shape suits it: it is
only ever consumed by the three methods paired with it.

[`SeriesVisualisation`](@ref) covers every curve-fit experiment and is the default.
"""
abstract type ResultVisualisation end

"""
    SeriesVisualisation()

Observed points with error bars plus a fitted line, one colour per group. The default,
and what relaxation, TRACT, nutation, diffusion and kinetics all use.
"""
struct SeriesVisualisation <: ResultVisualisation end

visualisationtype(::Experiment1D) = SeriesVisualisation()

# Trait forwarding: callers use the three-argument forms and never name a strategy.
function completeresultstate!(state, expt::Experiment1D)
    return completeresultstate!(state, expt, visualisationtype(expt))
end
resultpanel!(gui, state, expt::Experiment1D) = resultpanel!(gui, state, expt,
                                                            visualisationtype(expt))
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
fmt(x::Measurement, digits=4) = string(round(Measurements.value(x); sigdigits=digits),
                                       " ± ",
                                       round(Measurements.uncertainty(x); sigdigits=2))
fmt(x::Real, digits=4) = string(round(x; sigdigits=digits))
fmt(::Nothing, digits=4) = "n/a"

# Units and display names for parameters, fitted and derived alike. Derived quantities are
# *stored* in the unit named here (see `RegionResult`), so one number serves the summary
# and any tabular export. The 2D side keeps an equivalent table in `gui2d/summary.jl`
# (`PARAM_LABELS`); merging the two is the first step of the open question in PLAN.md about
# a canonical output format across 1D and 2D.
# Keys are unique across every experiment - a symbol means one quantity in one unit
# everywhere, which is why TRACT's cross-correlated rate is `:ηxy` and diffusion's
# viscosity is `:viscosity` rather than both being `:η`.
const PARAM_UNITS = Dict(:R => " s⁻¹",
                         :ν => " Hz",
                         :k => " s⁻¹",
                         :ηxy => " s⁻¹",
                         :τc => " ns",
                         :pulse90 => " µs",
                         :D => " ×10⁻¹⁰ m² s⁻¹",
                         :rH => " Å",
                         :viscosity => " mPa s",
                         :relative => " %")

# `:A` gets a global label - "Amplitude" is accurate for every fit here, whatever
# experiment it comes from. `:R` does not: relaxation and TRACT genuinely mean a
# relaxation rate by it, but nutation's damped-sinusoid model uses the same bare symbol
# for its decay rate, which "Relaxation rate" would mislabel. So `:R` stays out of this
# global table and is instead overridden per experiment, below, only where it is true.
const PARAM_LABELS = Dict(:A => "Amplitude",
                          :pulse90 => "90°",
                          :inhomogeneity => "B₁ inhom.",
                          :ηxy => "CCR rate (η)",
                          :τc => "Correlation time (τc)",
                          :STD_AF0 => "STD-AF₀",
                          :STD_AF => "STD-AF",
                          :STD_AF_max => "STD-AF_max",
                          :relative => "epitope",
                          :viscosity => "η")

"""
    paramlabel(expt, name) -> String
    paramunit(expt, name) -> String

Display name and unit for parameter `name` (fitted or derived). Experiment-dispatched,
not a single global table, precisely because a bare symbol can mean different things in
different experiments (see `PARAM_LABELS`'s note on `:R`) - the default falls back to the
symbols that genuinely are universal; an experiment overrides only the ones that aren't.
"""
paramlabel(::Experiment1D, name::Symbol) = get(PARAM_LABELS, name, string(name))
paramunit(::Experiment1D, name::Symbol) = get(PARAM_UNITS, name, "")

"""
    groupheader(expt, group) -> String

Display name for a grouping key (TRACT's `(which=:trosy,)`, STD's `(sat=:methyl,)`), used
as the bold header introducing that group's block in the results panel. Defaults to
[`groupname`](@ref)'s raw rendering of the values; TRACT overrides it to the same
TROSY/anti-TROSY wording its plot legend already uses (`seriesnames`).
"""
groupheader(::Experiment1D, group::NamedTuple) = groupname(group)

"""
    derivedheader(expt, r) -> String

Display name for the bold header introducing `r`'s derived-quantity block. Defaults to
[`groupheader`](@ref) when `r` is grouped, or "Derived" when it isn't - but a derived
quantity does not always belong to the group it happens to be recorded on (TRACT's η/τc
describe the *pair*, not specifically the TROSY series they are stored on for lack of a
region-level home - see `postfitglobal!` in `expt-tract.jl`), so this is a distinct hook
an experiment can override rather than reusing `groupheader` unconditionally.
"""
function derivedheader(e::Experiment1D, r::RegionResult)
    return isempty(r.group) ? "Derived" : groupheader(e, r.group)
end

"""
    paramblock(io, expt, params, width=nothing)

Write every entry of a parameter dictionary as `name value unit`, in insertion order, the
names padded to `width` (computed from `params` alone when omitted) so the values line
up. Callers building a whole panel from several dictionaries - several fitted-parameter
blocks, several derived blocks - pass a `width` computed once across all of them
([`panelwidth`](@ref)), so every block's values land in the same column; a lone caller
(`summary.txt`) leaves it to align to its own single block.
"""
function paramblock(io::IO, expt::Experiment1D, params, width=nothing)
    isempty(params) && return nothing
    w = something(width,
                  maximum(length(paramlabel(expt, name)) for name in keys(params)) + 2)
    for (name, value) in params
        println(io, "$(rpad(paramlabel(expt, name), w))$(fmt(value))$(paramunit(expt, name))")
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

The label-column width shared by every block in the results panel: the longest display
label among every fitted parameter and derived quantity shown for `activelabel`, plus a
gap, or `nothing` if there is nothing to show yet (leaving each block to size itself the
one time that matters least - before there is anything to align). Computed once so
"Amplitude" in the TROSY block and "Correlation time (τc)" in TRACT's results share one
column rather than each block aligning only to its own labels.
"""
function panelwidth(expt::Experiment1D, result, activelabel::AbstractString)
    len = 0
    for r in result
        r.region == activelabel || continue
        for k in keys(r.parameters)
            len = max(len, length(paramlabel(expt, k)))
        end
        r.postfitted || continue
        for k in keys(r.postparameters)
            len = max(len, length(paramlabel(expt, k)))
        end
    end
    return len == 0 ? nothing : len + 2
end

"""
    boldheader(text) -> RichText

A group/section heading in the results panel: bold, colon-suffixed, on its own line. The
one piece of visual hierarchy the panel uses - everything else in it is plain, aligned
parameter text.
"""
boldheader(text::AbstractString) = rich(rich(text * ":"; font=:bold), "\n")

"""
    plaintext(text) -> RichText

Wrap `text` in an explicit `font=:regular` span. `RichText`'s font is not scoped to each
sibling - a plain `String` child simply inherits whatever font the previous sibling left
active - so without this, the parameter block immediately after a `boldheader` would
render bold too, the bold state leaking straight past the header it belongs to.
"""
plaintext(text::AbstractString) = rich(text; font=:regular)

"""
    BLANK_RICHTEXT

Placeholder for "nothing to show". A genuinely empty `RichText` - `rich()`, zero children -
renders zero glyphs, and Makie's `GlyphCollection` cannot build itself from a zero-length
glyph vector (it can't infer the vector's `rotations` field as `Vector{Quaternionf}` from
an empty comprehension, and errors instead of rendering blank). A single space has exactly
one glyph, so it sidesteps that without being visible - not a cosmetic choice, without it
`secondarytext` crashes outright for every experiment that derives nothing (relaxation,
kinetics: `postfit!`'s default never marks a result `postfitted`, so its span list is
always empty), and `resultsheader` for a `NoFitting` series.
"""
const BLANK_RICHTEXT = rich(" ")

"""
    richtext(spans) -> RichText

`rich(spans...)`, substituting [`BLANK_RICHTEXT`](@ref) when `spans` is empty - see its
docstring for why that substitution is load-bearing, not decorative.
"""
richtext(spans) = isempty(spans) ? BLANK_RICHTEXT : rich(spans...)

"""
    resultsheader(expt, result, activelabel, width=nothing) -> RichText

The active region's raw fitted parameters (e.g. amplitude, rate) - the "primary" half of
[`summarytext`](@ref), split out so the GUI can show it and the derived
[`secondarytext`](@ref) (e.g. TRACT's τc) as visually separate blocks. Doesn't repeat the
region name, group headers only ([`groupheader`](@ref), e.g. TRACT's "TROSY"/"Anti-TROSY"),
since the region itself is already named elsewhere in the panel. A group header is bold,
matching [`secondarytext`](@ref)'s.

`width` aligns every block's values to the same column - see [`panelwidth`](@ref); omitted,
each block aligns only to its own labels.
"""
function resultsheader(e::Experiment1D, result, activelabel::AbstractString, width=nothing)
    spans = Any[]
    for r in result
        r.region == activelabel || continue
        block = paramtext(e, r.parameters, width)
        isempty(block) && continue
        isempty(r.group) || push!(spans, boldheader(groupheader(e, r.group)))
        push!(spans, plaintext(block), "\n")
    end
    return richtext(spans)
end

"""
    secondarytext(expt, result, activelabel, width=nothing) -> RichText

The active region's derived quantities alone (TRACT's τc and η, nutation's 90° pulse,
diffusion's D and rH) - empty when the experiment derives none, or when the active region
isn't ready yet (e.g. only one of a TRACT pair fitted so far). Like `resultsheader`, it
just reads `RegionResult.postparameters`, whatever wrote them - but its header
([`derivedheader`](@ref)) is not simply `resultsheader`'s, since a derived quantity does
not always belong to the group it happens to be recorded on.

`width` as in `resultsheader`.
"""
function secondarytext(e::Experiment1D, result, activelabel::AbstractString, width=nothing)
    spans = Any[]
    for r in result
        (r.region == activelabel && r.postfitted) || continue
        block = paramtext(e, r.postparameters, width)
        isempty(block) && continue
        push!(spans, boldheader(derivedheader(e, r)))
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
        header = isempty(r.group) ? r.region : "$(r.region) ($(groupheader(e, r.group)))"
        println(io, header)
        println(io, "-"^length(header))
        paramblock(io, e, r.parameters)
        r.postfitted && paramblock(io, e, r.postparameters)
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
    # Most experiments label each curve directly (group counts/names vary, e.g. STD's
    # saturation frequencies); experiments with a fixed, small set of named groups (e.g.
    # TRACT's TROSY/anti-TROSY) get a proper axislegend instead via `seriesnames`.
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
