# Makie plotting shared by every experiment type: GLMakie on screen, CairoMakie for the
# saved PDFs (see `_save_results`). A data panel above a residual panel a third its height,
# sharing an x-axis, drawn in the default Makie palette.

"Default Makie palette."
const PALETTE = Makie.wong_colors()
palettecolor(i) = PALETTE[mod1(i, length(PALETTE))]

"Axis attributes for a plot tiled with many others by [`combineplots`](@ref)."
const COMPACT = (titlesize=8, xlabelsize=7, ylabelsize=7, xticklabelsize=6,
                 yticklabelsize=6)

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
    for (i, expt) in enumerate(experiments)
        row, col = fldmod1(i, ncols)
        gl = f[row, col] = GridLayout()
        plotresult!(gl, expt, result; axiskw=COMPACT)
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
