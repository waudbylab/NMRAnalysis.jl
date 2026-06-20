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

Collect `param` (default: the experiment's [`primaryparam`](@ref)) for every peak
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

    pcol = findfirst(==(string(param)), header)
    isnothing(pcol) &&
        error("Parameter \"$param\" not found in $file. Available: $(available_params(path))")
    ecol = findfirst(==("$(param)_err"), header)
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

function _resolveparam(source, param)
    source isa AbstractVector && return _resolveparam(first(source), param)
    param isa Symbol && return param
    return _defaultparam(source)
end

_todatasets(source::AbstractVector, param; kw...) =
    [_onedataset(s, param; kw...) for s in source]
_todatasets(source, param; kw...) = [_onedataset(source, param; kw...)]

_onedataset(s::FixedPeakExperiment, param; kw...) = summary_dataset(s, param; kw...)
_onedataset(s::AbstractString, param; kw...) = summary_dataset(s, param; kw...)

# --- plotting --------------------------------------------------------------

"""
    summaryplot(source, param=<default>; ylabel, title, include_unassigned)

Plot a fitted parameter against residue number.

`source` may be a live experiment, a saved `results.csv` (or its folder), or a
vector of any of these. `param` defaults to the experiment's
[`primaryparam`](@ref) (or the first derived column of a file).

- With backbone/amide labels → a scatter plot of value vs residue number, with
  error bars.
- When any atom-typed labels are present (e.g. methyls `I13CD1`, `L26CD2`) → a
  bar plot ordered by `(residue, atom)` with peak labels as tick marks, so
  stereospecific pairs don't overlap.
- A vector `source` → vertically stacked panels sharing the residue axis.
- Unassigned peaks (the default `X#` names) are omitted unless every peak is
  unassigned, or `include_unassigned=true`.

Uses whichever Makie backend is active and returns the `Figure`, so the result
displays interactively under GLMakie and can be saved with
`save("summary.pdf", fig)` under CairoMakie.
"""
function summaryplot(source, param=nothing; ylabel=nothing, title="",
                     include_unassigned=false)
    param = _resolveparam(source, param)
    datasets = _todatasets(source, param; include_unassigned=include_unassigned)

    # consistent plot kind across all panels
    usebar = any(ds -> has_atom_labels(p.label for p in ds.points), datasets)
    ylab = something(ylabel, paramlabel(param))

    fig = Figure()
    axes = Axis[]
    n = length(datasets)
    for (i, ds) in enumerate(datasets)
        ax = Axis(fig[i, 1];
                  xlabel=(i == n ? (usebar ? "" : "Residue number") : ""),
                  ylabel=ylab,
                  title=(n > 1 ? ds.name : title))
        _drawdataset!(ax, ds, usebar)
        push!(axes, ax)
    end
    n > 1 && !usebar && linkxaxes!(axes...)

    return fig
end

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
        plt = barplot!(ax, xs, ys; color=(:steelblue, 0.8))
        errorbars!(ax, xs, ys, es; whiskerwidth=6)
        ax.xticks = (xs, [p.label for p in points])
        ax.xticklabelrotation = π / 2
        ax.xticklabelsize = 9
        return plt
    else
        xs = [Float64(p.resnum) for p in points]
        ys = [p.value for p in points]
        es = [p.uncertainty for p in points]
        errorbars!(ax, xs, ys, es; whiskerwidth=6, color=:grey40)
        plt = scatter!(ax, xs, ys; color=:steelblue)
        return plt
    end
end
