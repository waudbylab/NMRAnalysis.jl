"""
    experimentgroups(key, experiments) -> Vector{Vector}

Group `experiments` by the value of `key(experiment)`, preserving the order
in which each group's key was first encountered. Used by `overlayplots` to
find sets of experiments that share a saturation power or saturation time
(and therefore differ only in the other one), so they can be drawn together.
"""
function experimentgroups(key, experiments)
    order = Any[]
    groups = Dict{Any,Vector{eltype(experiments)}}()
    for expt in experiments
        k = key(expt)
        if !haskey(groups, k)
            groups[k] = eltype(experiments)[]
            push!(order, k)
        end
        push!(groups[k], expt)
    end
    return [groups[k] for k in order]
end

"""
Solid, lightened marker fill colours for the *data* panel of an overlay
plot, one per curve (indexed the same way as the full-strength `color=i`
used for fit lines and for the residual panel) — not transparency: alpha
fades the error bars along with the marker fill, and overlapping
semi-transparent error bars from several experiments composite into a muddy
mess. A flat lighter colour keeps error bars crisp (their stroke colour is
still forced to black — see the `markerstrokecolor` on each data scatter!)
while letting the full-strength fit line read as the most prominent element
on the plot. The residual panel has no competing fit line, so its markers
use the regular full-strength `color=i` instead.
"""
const _OVERLAY_MARKER_COLORS = [:lightskyblue, :navajowhite, :palegreen,
                                :lightpink, :plum, :lightgray]
overlaymarkercolor(i) = _OVERLAY_MARKER_COLORS[mod1(i, length(_OVERLAY_MARKER_COLORS))]

"""
    overlayplots(result::FitResult) -> Vector

Additional plots overlaying similar experiments on shared axes, so that
trends across a series are easy to compare directly (issue #39):

- all off-resonance R1ρ experiments together, one curve per spin-lock power
- CEST experiments sharing a saturation power (ν1), one curve per
  saturation time
- CEST experiments sharing a saturation time, one curve per saturation
  power

These are kept separate from `combineplots(plot(result))` — the per-experiment
grid — rather than folded into it, so that grid doesn't grow denser (and
illegible) as overlays are added; save or display them as their own figure(s)
instead. Returns an empty vector if there is nothing to usefully overlay
(fewer than two experiments in a group).
"""
function overlayplots(result::FitResult)
    prob = result.prob
    experiments = prob.experiments
    plots = Any[]

    offres = filter(e -> e isa R1rhoOffResExperiment, experiments)
    cest = filter(e -> e isa CESTExperiment, experiments)
    (length(offres) > 1 || length(cest) > 1) || return plots

    simulate!(prob, result.params_value)  # ensure predicted_intensities are current
    params_value = result.params_value

    length(offres) > 1 && push!(plots, overlayr1rhooffres(offres, params_value))

    if length(cest) > 1
        for group in experimentgroups(e -> round(e.ν1; digits=1), cest)
            length(group) > 1 &&
                push!(plots, overlaycest(group, :saturation_time, params_value))
        end
        for group in experimentgroups(e -> round(e.saturation_time; digits=6), cest)
            length(group) > 1 && push!(plots, overlaycest(group, :ν1, params_value))
        end
    end

    return plots
end

"""
    overlayr1rhooffres(experiments, params_value) -> Plot

Overlay all off-resonance R1ρ experiments on shared axes, one colour per
spin-lock power, sorted by offset. Data, fit, and residuals for a given
experiment share a colour; residual bands are shaded grey (rather than the
green used for a single experiment's plot) so they read clearly with
multiple overlapping colours. The fitted chemical shift(s)
(`params_value.spin.delta`) are marked with dashed vertical lines, as on
the individual per-experiment plot.

The legend is drawn inset (`:topright`), not outside the axes: an outer
legend only widens the data panel, not the residual panel below it, which
throws off their alignment under the shared `link=:x`.

Draw order (bottom to top): resonance-position vlines, then residual bands
and zero-line, then data (markers + error bars) for every experiment, then
every fit line in a second pass — so the fit lines render on top of every
experiment's data rather than only their own, and the data markers render
above the residual bands. See `_OVERLAY_MARKER_COLORS`.
"""
function overlayr1rhooffres(experiments, params_value)
    delta = params_value.spin.delta
    sorted = sort(experiments; by=expt -> expt.νSL)

    p1 = plot(; frame=:box, xlabel="Spin-lock offset / ppm", ylabel="R₁ρ / s⁻¹",
              title="Off-resonance R₁ρ (combined)", grid=nothing, xflip=true,
              legend=:topright)
    p2 = plot(; frame=:box, xlabel="Spin-lock offset / ppm", ylabel="Residual / σ",
              grid=nothing, xflip=true, legend=nothing)

    hspan!(p2, [-2, 2]; color=:grey, alpha=0.3, lw=0, la=0, primary=false)
    hspan!(p2, [-1, 1]; color=:grey, alpha=0.5, lw=0, la=0, primary=false)

    hline!(p2, [0]; color=:black, lw=0.5, primary=false)

    maxwres = 0.0
    curves = map(enumerate(sorted)) do (i, expt)
        x = expt.offsets_ppm
        sortidx = sortperm(x)
        yobs = expt.observed_intensities
        ypred = expt.predicted_intensities
        wres = residuals(expt)
        maxwres = max(maxwres, maximum(abs, wres))
        label = "$(Int(round(expt.νSL; digits=0))) Hz"
        return (; i, x, sortidx, yobs, ypred, wres, label)
    end

    for c in curves
        scatter!(p1, c.x, c.yobs; ms=3, markercolor=overlaymarkercolor(c.i),
                 markerstrokecolor=:black, label=c.label)
        scatter!(p2, c.x, c.wres; ms=3, color=c.i, label=nothing)
    end
    for c in curves
        plot!(p1, c.x[c.sortidx], c.ypred[c.sortidx]; lw=1.5, color=c.i, label=nothing)
    end

    hline!(p1, [0.0]; color=:black, lw=0.5, primary=false, z_order=:back)
    vline!(p1, delta; ls=:dash, color=:black, label=nothing, z_order=:back)
    vline!(p2, delta; ls=:dash, color=:black, label=nothing, z_order=:back)

    ylims!(p2, -maxwres * 1.2, maxwres * 1.2)

    return plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]), link=:x)
end

"""
    overlaycest(experiments, varyby::Symbol, params_value) -> Plot

Overlay CEST experiments that share every acquisition parameter except
`varyby` (`:saturation_time` or `:ν1`) on shared axes, sorted by `varyby`,
one colour per curve. Data, fit, and residuals for a given experiment share
a colour; residual bands are shaded grey rather than green, as in
`overlayr1rhooffres`. The fitted chemical shift(s) (`params_value.spin.delta`)
are marked with dashed vertical lines, as on the individual per-experiment
plot; observed/predicted intensities are normalised by the fitted `I0` (see
`CESTExperiment`'s `plot_result`) so the baseline sits at 1.0.

The legend is drawn inset (`:bottomright`, where the CEST dip is not) rather
than outside the axes, for the same alignment reason as
`overlayr1rhooffres`.

Draw order (bottom to top): resonance-position vlines, then residual bands
and zero-line, then data (markers + error bars) for every experiment, then
every fit line in a second pass. See `_OVERLAY_MARKER_COLORS`.
"""
function overlaycest(experiments, varyby::Symbol, params_value)
    delta = params_value.spin.delta
    sorted = sort(experiments; by=expt -> getproperty(expt, varyby))
    fixed = sorted[1]
    title = if varyby == :saturation_time
        "CEST ($(Int(round(fixed.ν1; digits=0))) Hz, combined)"
    else
        "CEST ($(Int(round(fixed.saturation_time * 1000; digits=0))) ms, combined)"
    end

    p1 = plot(; frame=:box, ylabel="Normalised intensity", title=title, grid=nothing,
              xflip=true, legend=:bottomright)
    p2 = plot(; frame=:box, xlabel="Saturation frequency (ppm)", ylabel="Residual / σ",
              grid=nothing, xflip=true, legend=nothing)

    hspan!(p2, [-2, 2]; color=:grey, alpha=0.3, lw=0, la=0, primary=false)
    hspan!(p2, [-1, 1]; color=:grey, alpha=0.5, lw=0, la=0, primary=false)
    hline!(p2, [0]; color=:black, lw=0.5, primary=false)

    maxwres = 0.0
    curves = map(enumerate(sorted)) do (i, expt)
        x = expt.δsat
        sortidx = sortperm(x)
        I0 = params_value.nuisance[Symbol("CEST_", field_label(expt), "_I0")]
        yobs = expt.observed_intensities ./ I0
        ypred = expt.predicted_intensities ./ I0
        wres = residuals(expt)
        maxwres = max(maxwres, maximum(abs, wres))
        label = varyby == :saturation_time ?
                "$(Int(round(expt.saturation_time * 1000; digits=0))) ms" :
                "$(Int(round(expt.ν1; digits=0))) Hz"
        return (; i, x, sortidx, yobs, ypred, wres, label)
    end

    for c in curves
        scatter!(p1, c.x[c.sortidx], c.yobs[c.sortidx]; ms=3,
                 markercolor=overlaymarkercolor(c.i), markerstrokecolor=:black,
                 label=c.label)
        scatter!(p2, c.x[c.sortidx], c.wres[c.sortidx]; ms=3, color=c.i, label=nothing)
    end
    for c in curves
        plot!(p1, c.x[c.sortidx], c.ypred[c.sortidx]; lw=1.5, color=c.i, label=nothing)
    end

    hline!(p1, [0.0]; color=:black, lw=0.5, primary=false, z_order=:back)
    vline!(p1, delta; ls=:dash, color=:black, label=nothing, z_order=:back)
    vline!(p2, delta; ls=:dash, color=:black, label=nothing, z_order=:back)

    ylims!(p2, -maxwres * 1.2, maxwres * 1.2)

    return plot(p1, p2; layout=grid(2, 1; heights=[0.75, 0.25]), link=:x)
end
