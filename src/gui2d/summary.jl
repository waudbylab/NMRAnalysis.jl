# Results-summary plotting: parameter vs residue number.
#
# `summaryplot` is a plain Makie plot built with whichever backend is currently
# active (interactive under GLMakie, saveable under CairoMakie). It takes a live
# experiment, a saved `results.csv` (or its folder), or a vector of these, and
# returns a standard Makie figure the user can `save`/display as usual.

"""
    SummaryPoint

One peak's contribution to a summary plot: identity plus a single value and its
uncertainty.
"""
struct SummaryPoint
    label::String
    resnum::Int
    resname::Char
    atom::String
    value::Float64
    uncertainty::Float64
end

"""
    SummaryDataset

A named collection of [`SummaryPoint`](@ref)s for one parameter, from one source
(experiment or file).
"""
struct SummaryDataset
    name::String
    param::Symbol
    points::Vector{SummaryPoint}
end

# Default y-axis labels for known parameters (deliberately makes no R1/R2
# assumption for the generic exponential rate :R).
const PARAM_LABELS = Dict(:hetnoe => "Heteronuclear NOE",
                          :R20 => "R₂⁰ / s⁻¹",
                          :PRE => "Γ₂ / s⁻¹",
                          :eta => "η / s⁻¹",
                          :R => "Relaxation rate / s⁻¹",
                          :R1 => "R₁ / s⁻¹",
                          :R2 => "R₂ / s⁻¹")

paramlabel(param) = get(PARAM_LABELS, Symbol(param), string(param))

_resultsfile(path) = isdir(path) ? joinpath(path, "results.csv") : path

function _sourcename(path::AbstractString)
    b = basename(isdir(path) ? path : dirname(path))
    return isempty(b) ? string(path) : b
end
_sourcename(expt) = string(nameof(typeof(expt)))

# --- value extraction ------------------------------------------------------

"Value/uncertainty for `param` from a live peak (post-fit params first, then raw)."
function _peakvalue(peak, param::Symbol)
    if haskey(peak.postparameters, param)
        p = peak.postparameters[param]
    elseif haskey(peak.parameters, param)
        p = peak.parameters[param]
    else
        return nothing
    end
    return (Float64(to_value(p.value[][1])), Float64(to_value(p.uncertainty[][1])))
end

"Drop unassigned peaks (resnum ≤ 0), unless that would leave nothing to plot."
function _filterunassigned(points, include_unassigned)
    include_unassigned && return points
    assigned = filter(p -> p.resnum > 0, points)
    return isempty(assigned) ? points : assigned
end

# --- dataset construction --------------------------------------------------

"""
    summary_dataset(expt, param; name, include_unassigned) -> SummaryDataset

Collect `param` (default: the experiment's `primaryparam`) for every peak
in a live experiment.
"""
function summary_dataset(expt::FixedPeakExperiment, param::Symbol=primaryparam(expt);
                         name="", include_unassigned=false)
    points = SummaryPoint[]
    for peak in expt.peaks[]
        vu = _peakvalue(peak, param)
        isnothing(vu) && continue
        lbl = parse_label(peak.label[])
        push!(points,
              SummaryPoint(peak.label[], lbl.resnum, lbl.onelettercode, lbl.atom,
                           vu[1], vu[2]))
    end
    points = _filterunassigned(points, include_unassigned)
    return SummaryDataset(isempty(name) ? _sourcename(expt) : name, param, points)
end

"""
    summary_dataset(path, param; name, include_unassigned) -> SummaryDataset

Read `param` from a saved `results.csv` (or a folder containing one). The label,
residue number and atom are re-derived from each row's label.
"""
function summary_dataset(path::AbstractString, param::Symbol;
                         name="", include_unassigned=false)
    file = _resultsfile(path)
    header, rows = _readtable(file)

    # Arrayed parameters are written per plane as e.g. amp[1], amp[2]; if the
    # bare name (`amp`) is requested, default to the first index (`amp[1]`),
    # matching how a live experiment plots the first slice.
    pcol = findfirst(==(string(param)), header)
    isnothing(pcol) && (pcol = findfirst(==("$(param)[1]"), header))
    isnothing(pcol) &&
        error("Parameter \"$param\" not found in $file. Available: $(available_params(path))")
    pname = header[pcol]
    ecol = findfirst(==("$(pname)_err"), header)
    lcol = something(findfirst(==("label"), header), 1)

    points = SummaryPoint[]
    for row in rows
        length(row) < pcol && continue
        v = tryparse(Float64, row[pcol])
        isnothing(v) && continue
        u = (isnothing(ecol) || ecol > length(row)) ? 0.0 :
            something(tryparse(Float64, row[ecol]), 0.0)
        label = row[lcol]
        lbl = parse_label(label)
        push!(points, SummaryPoint(label, lbl.resnum, lbl.onelettercode, lbl.atom, v, u))
    end
    points = _filterunassigned(points, include_unassigned)
    return SummaryDataset(isempty(name) ? _sourcename(path) : name, param, points)
end

"""
    _readtable(file) -> (header, rows)

Read a `results.csv`: skip blank and `#`-comment lines, take the first remaining
line as the (comma- or whitespace-separated) header, the rest as rows of strings.
"""
function _readtable(file)
    header = String[]
    rows = Vector{String}[]
    for line in eachline(file)
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        fields = String.(splitfields(s))
        if isempty(header)
            header = fields
        else
            push!(rows, fields)
        end
    end
    return header, rows
end

# --- available parameters / defaults ---------------------------------------

"""
    available_params(source) -> Vector{Symbol}

Parameters that can be plotted from a live experiment (post-fit and raw
parameter keys) or a saved file (value columns, excluding identity columns).
"""
function available_params(expt::FixedPeakExperiment)
    isempty(expt.peaks[]) && return Symbol[]
    p = first(expt.peaks[])
    return unique([collect(keys(p.postparameters)); collect(keys(p.parameters))])
end

const _IDENTITY_COLS = ("label", "resnum", "resname", "atom")

function available_params(path::AbstractString)
    header, _ = _readtable(_resultsfile(path))
    return [Symbol(h) for h in header if !endswith(h, "_err") && !(h in _IDENTITY_COLS)]
end

_defaultparam(expt::FixedPeakExperiment) = primaryparam(expt)

function _defaultparam(path::AbstractString)
    header, _ = _readtable(_resultsfile(path))
    fixed = ("label", "resnum", "resname", "atom", "x", "y", "R2x", "R2y")
    for h in header
        (endswith(h, "_err") || h in fixed || startswith(h, "amp[")) && continue
        return Symbol(h)
    end
    return Symbol("amp[1]")
end

_onedataset(s::FixedPeakExperiment, param; kw...) = summary_dataset(s, param; kw...)
_onedataset(s::AbstractString, param; kw...) = summary_dataset(s, param; kw...)

# Resolve one parameter per source. `param` may be:
#   nothing                  -> each source's own default (primaryparam / first
#                               derived column), so heterogeneous experiments
#                               (e.g. relaxation + hetNOE) each plot their own;
#   a Symbol                 -> the same parameter for every source;
#   a vector (Symbol/nothing)-> one per source, nothing meaning that source's default.
function _paramper(sources, param)
    if param === nothing
        return [_defaultparam(s) for s in sources]
    elseif param isa Symbol
        return fill(param, length(sources))
    elseif param isa AbstractVector
        length(param) == length(sources) ||
            error("Got $(length(param)) parameters for $(length(sources)) sources")
        return [p === nothing ? _defaultparam(s) : Symbol(p) for (s, p) in zip(sources, param)]
    else
        error("`param` must be a Symbol, a vector of Symbols, or nothing")
    end
end

function _ylabelper(params, ylabel)
    if ylabel === nothing
        return [paramlabel(p) for p in params]
    elseif ylabel isa AbstractVector
        length(ylabel) == length(params) ||
            error("Got $(length(ylabel)) ylabels for $(length(params)) panels")
        return [y === nothing ? paramlabel(p) : string(y) for (p, y) in zip(params, ylabel)]
    else
        return fill(string(ylabel), length(params))
    end
end

# --- plotting --------------------------------------------------------------

"""
    summaryplot(source; param=<default>, ylabel, title, size, include_unassigned)
    summaryplot(source1, source2, ...; kwargs...)

Plot a fitted parameter against residue number.

`source` may be a live experiment, a saved `results.csv` (or its folder), or a
vector of any of these (which gives vertically stacked panels). Multiple sources
may also be passed as separate positional arguments instead of a vector.

`param` selects which parameter to plot:
- omitted/`nothing` → each source's own default (its `primaryparam` or
  first derived column), so a mix of experiment types — e.g. relaxation and
  hetNOE — each plot their own result;
- a single `Symbol` → the same parameter for every source;
- a vector → one parameter per source (a `nothing` entry uses that source's
  default).

`ylabel` likewise may be a single label applied to all panels or a vector of
per-panel labels; by default each panel is labelled from its parameter.

`size` sets the figure size in pixels, e.g. `size=(800, 400)`.

- Backbone/amide labels → a scatter of value vs residue number with error bars.
- Any atom-typed labels (e.g. methyls `I13CD1`, `L26CD2`) → a bar plot ordered
  by `(residue, atom)` with peak labels as ticks, so stereospecific pairs don't
  overlap. Decided per panel.
- Unassigned peaks (the default `X#` names) are omitted unless every peak is
  unassigned, or `include_unassigned=true`.

Uses whichever Makie backend is active and returns the `Figure`, so the result
displays interactively under GLMakie and can be saved with
`save("summary.pdf", fig)` under CairoMakie.
"""
function summaryplot(source; param=nothing, ylabel=nothing, title="", size=nothing,
                     include_unassigned=false)
    sources = source isa AbstractVector ? collect(source) : [source]
    params = _paramper(sources, param)
    ylabels = _ylabelper(params, ylabel)

    datasets = [_onedataset(s, p; include_unassigned=include_unassigned)
                for (s, p) in zip(sources, params)]

    figkw = isnothing(size) ? NamedTuple() : (; size=size)
    fig = Figure(; figkw...)
    axes = Axis[]
    n = length(datasets)
    panelbar = [has_atom_labels(p.label for p in ds.points) for ds in datasets]
    for (i, ds) in enumerate(datasets)
        usebar = panelbar[i]
        ax = Axis(fig[i, 1];
                  xlabel=(i == n ? (usebar ? "" : "Residue number") : ""),
                  ylabel=ylabels[i],
                  title=(n > 1 ? ds.name : title),
                  xgridvisible=false, ygridvisible=false)
        _drawdataset!(ax, ds, usebar)
        hlines!(ax, [0]; linewidth=0)  # invisible: forces zero into the y-range
        push!(axes, ax)
    end
    # share the residue axis across stacked scatter panels
    n > 1 && !any(panelbar) && linkxaxes!(axes...)

    return fig
end

# Convenience varargs form: summaryplot("a/", "b/", "c/"; kwargs...)
summaryplot(s1, s2, rest...; kw...) = summaryplot([s1, s2, rest...]; kw...)

function _drawdataset!(ax, ds::SummaryDataset, usebar)
    points = ds.points
    if isempty(points)
        return scatter!(ax, Float64[], Float64[])
    end
    if usebar
        order = sortperm([(p.resnum, p.atom) for p in points])
        points = points[order]
        xs = collect(1:length(points))
        ys = [p.value for p in points]
        es = [p.uncertainty for p in points]
        plt = barplot!(ax, xs, ys)
        errorbars!(ax, xs, ys, es; whiskerwidth=6, color=:black)
        ax.xticks = (xs, [p.label for p in points])
        ax.xticklabelrotation = π / 4
        return plt
    else
        xs = [Float64(p.resnum) for p in points]
        ys = [p.value for p in points]
        es = [p.uncertainty for p in points]
        errorbars!(ax, xs, ys, es; whiskerwidth=6, color=:black)
        plt = scatter!(ax, xs, ys; color=:steelblue)
        return plt
    end
end
