# Makie plotting shared by every experiment type: GLMakie on screen, CairoMakie for the
# saved PDFs (see `_save_results`). A data panel above a residual panel a third its height,
# sharing an x-axis, drawn in the default Makie palette.

"Default Makie palette."
const PALETTE = Makie.wong_colors()
palettecolor(i) = PALETTE[mod1(i, length(PALETTE))]

"Default Makie `Axis` font sizes, used by [`combineplots`](@ref) while there are few
enough panels for the default sizing to stay legible."
const DEFAULTFONTS = (titlesize=16, xlabelsize=16, ylabelsize=16, xticklabelsize=16,
                      yticklabelsize=16)

"The smallest [`combineplots`](@ref) ever shrinks its fonts to, for a grid of many panels."
const COMPACTFONTS = (titlesize=12, xlabelsize=10.5, ylabelsize=10.5, xticklabelsize=9,
                      yticklabelsize=9)

"""
    combinedfontscale(n) -> NamedTuple

Axis font sizes for a [`combineplots`](@ref) grid of `n` panels: [`DEFAULTFONTS`](@ref) up
to a single row, shrinking down to [`COMPACTFONTS`](@ref) as the grid grows to several rows.
"""
function combinedfontscale(n)
    t = clamp((n - 4) / (20 - 4), 0, 1)
    return map((lo, hi) -> lo + t * (hi - lo), DEFAULTFONTS, COMPACTFONTS)
end

"""
    resultpanels!(gl; xlabel, ylabel, title="", xreversed=false, axiskw=(;))

The two panels of a result plot: data above, weighted residuals below at a third the
height, sharing an x-axis. Returns both axes.
"""
function resultpanels!(gl; xlabel, ylabel, title="", xreversed=false, axiskw=(;))
    common = (; xgridvisible=false, ygridvisible=false, xreversed=xreversed, axiskw...)
    # Makie.Axis: ComponentArrays also exports `Axis`, so the bare name is ambiguous.
    # The upper panel takes no xlabel: the shared axis is labelled once, underneath.
    ax1 = Makie.Axis(gl[1, 1]; ylabel=ylabel, title=title, common...)
    ax2 = Makie.Axis(gl[2, 1]; xlabel=xlabel, ylabel="Residual / σ", common...)
    linkxaxes!(ax1, ax2)
    rowsize!(gl, 1, Auto(false, 3))
    rowsize!(gl, 2, Auto(false, 1))
    rowgap!(gl, 4)
    return ax1, ax2
end

"""
    residualbands!(ax; color=:limegreen)

Shade ±1σ and ±2σ on a residual panel and mark zero. Overlay plots pass `color=:grey`,
several curves of different colours reading badly against green.
"""
function residualbands!(ax; color=:limegreen)
    hspan!(ax, -2, 2; color=(color, 0.3))
    hspan!(ax, -1, 1; color=(color, 0.5))
    hlines!(ax, [0]; color=:black, linewidth=0.5)
    return ax
end

"""
    measured!(ax, x, y; color, label=nothing, strokecolor=nothing)

Plot measurements carrying uncertainties: an error bar and a marker per point. Makie draws
neither from a `Measurement` on its own.
"""
function measured!(ax, x, y; color, label=nothing, strokecolor=nothing)
    value = Measurements.value.(y)
    err = Measurements.uncertainty.(y)
    errorbars!(ax, x, value, err; color=:black, linewidth=0.7, whiskerwidth=4)
    return if isnothing(label)
        scatter!(ax, x, value; color=color, markersize=6,
                 strokecolor=something(strokecolor, color), strokewidth=0.5)
    else
        scatter!(ax, x, value; color=color, markersize=6,
                 strokecolor=something(strokecolor, color), strokewidth=0.5, label=label)
    end
end

"""
    resultfigure(expt, fitresult) -> Figure

One experiment's result plot as a figure of its own, for `experiments/<name>.pdf`.
"""
function resultfigure(expt, fitresult)
    f = Figure(; size=(500, 400))
    gl = f[1, 1] = GridLayout()
    plotresult!(gl, expt, fitresult)
    return f
end

"""
    combineplots(result) -> Figure

Every experiment of a fit on one grid, at most four columns wide, with the fonts and the
figure scaled so that a fit of many experiments stays legible.
"""
function combineplots(result::FitResult)
    experiments = result.prob.experiments
    simulate!(result.prob, result.params_value)
    n = length(experiments)
    ncols = clamp(n, 1, 4)
    nrows = ceil(Int, n / ncols)

    f = Figure(; size=(max(1200, ncols * 350), max(800, nrows * 280)))
    axiskw = combinedfontscale(n)
    for (i, expt) in enumerate(experiments)
        row, col = fldmod1(i, ncols)
        gl = f[row, col] = GridLayout()
        plotresult!(gl, expt, result; axiskw=axiskw)
    end
    return f
end

"""
    plotresults(result) -> Vector{Figure}

One figure per experiment, in the problem's own order.
"""
function plotresults(result::FitResult)
    simulate!(result.prob, result.params_value)
    return [resultfigure(expt, result) for expt in result.prob.experiments]
end

"""
    correlationheatmap(result::FitResult) -> Figure

Heatmap of the correlation matrix of the fitted (non-fixed) parameters, axes labelled with
the same display names as the parameter tables — the picture version of `correlation.csv`,
for spotting at a glance which parameters trade off against each other.
"""
function correlationheatmap(result::FitResult)
    state_labels = states(result.prob.model)
    fields = _unique_fields(result.prob.experiments)
    freeitems = _flatten_params_items(result.params)[result.freeidx]
    labels = [_pretty_label(item, state_labels, fields) for item in freeitems]
    n = length(labels)

    f = Figure(; size=(max(500, 55 * n + 220), max(450, 55 * n + 120)))
    ax = Makie.Axis(f[1, 1]; title="Parameter correlation",
                    xticks=(1:n, labels), yticks=(1:n, labels),
                    xticklabelrotation=π / 4, yreversed=true)
    hm = heatmap!(ax, 1:n, 1:n, result.cor; colormap=:RdBu, colorrange=(-1, 1))
    Colorbar(f[1, 2], hm; label="r")
    return f
end
