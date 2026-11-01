# Figures for saving, drawn at a fixed size in points whatever the window's size.

const PALETTE = Makie.wong_colors()

"The colour of multiplet `m`: black for the reference, then the Wong palette in order."
function multipletcolour(m, α=1)
    m == REF && return Makie.RGBAf(0, 0, 0, α)
    c = PALETTE[mod1(m - 1, length(PALETTE))]
    return Makie.RGBAf(c.r, c.g, c.b, α)
end

"Axis label for a nucleus such as \"1H\" or \"19F\", with the mass number as a superscript."
function shiftlabel(nucleus::AbstractString)
    m = match(r"^(\d+)(\D+)$", nucleus)
    isnothing(m) && return "$nucleus chemical shift (ppm)"
    return rich(Makie.superscript(String(m[1])), String(m[2]), " chemical shift (ppm)")
end

"""
    stackplot(p::PeakTracker; size=(252, 340), fontsize=8, markers=true) -> Figure

A publication figure of the chemical shift range last shown in the window, with its
vertical scale and any alignment. The size is in points, 8.9 × 12 cm by default, a single
journal column. There is no legend, so every figure has the same plot area; picks are
coloured by multiplet as in the window. Save it with CairoMakie for vector output:

```julia
save("stack.pdf", stackplot(p); backend=CairoMakie, pt_per_unit=1)
```
"""
function stackplot(p::PeakTracker; size=(252, 340), fontsize=8, markers=true)
    s, v, t = p.series, p.view, tracks(p)
    lims = p.axis.finallimits[]
    lo, hi = lims.origin[1], lims.origin[1] + lims.widths[1]
    offsets = v.offsets[]

    fig = Figure(; size, fontsize, figure_padding=(4, 10, 4, 4))
    ax = Axis(fig[1, 1]; xlabel=shiftlabel(s.nucleus),
              yticksvisible=false, yticklabelsvisible=false,
              xgridvisible=false, ygridvisible=false,
              leftspinevisible=false, rightspinevisible=false, topspinevisible=false,
              spinewidth=0.75, xtickwidth=0.75, xticksize=3)

    top = (length(s) - 1) * v.spacing
    for i in eachindex(s)
        x = shifts(s) .+ offsets[i]
        sel = findall(d -> lo ≤ d ≤ hi, x)
        isempty(sel) && continue
        y = ypos.(v.spacing, v.gain[], v.norm, i, intensities(s, i)[sel])
        top = max(top, maximum(y))
        lines!(ax, x[sel], y; color=:black, linewidth=0.5)
    end

    if markers
        for (l, m) in ((l, m) for m in eachindex(multiplets(t)) for l in linesof(t, m))
            qs = filter(q -> lo ≤ shown(v, q.plane, q.shift) ≤ hi, picksof(t, l))
            isempty(qs) && continue
            segments = Point2f[]
            for (a, c) in zip(qs, Iterators.drop(qs, 1))
                c.plane == a.plane + 1 && push!(segments, location(v, a), location(v, c))
            end
            isempty(segments) ||
                linesegments!(ax, segments; color=multipletcolour(m), linewidth=0.75)
            scatter!(ax, [location(v, q) for q in qs]; color=multipletcolour(m),
                     marker=:dtriangle, markersize=4)
        end
    end

    # limits given high-to-low: Makie takes the axis direction from their order, so this
    # is what makes chemical shift decrease to the right
    limits!(ax, hi, lo, -0.5v.spacing, top + 0.3v.spacing)
    return fig
end

"""
    multipletplot(p, m) -> Figure

Multiplet `m` against the first numeric coordinate, or against plane number where there is
none: each line's shift as an open marker, and the centre filled.
"""
function multipletplot(p::PeakTracker, m; size=(340, 250), fontsize=10)
    t, s = tracks(p), p.series
    name = findfirst(v -> eltype(v) <: Real, p.coordinates)
    x = isnothing(name) ? collect(1:length(s)) : Float64.(p.coordinates[name])
    fig = Figure(; size, fontsize)
    ax = Axis(fig[1, 1]; xlabel=isnothing(name) ? "plane" : string(name),
              ylabel=shiftlabel(s.nucleus), title=multiplets(t)[m].label)
    ls = orderedlines(t, m)
    if length(ls) > 1
        for l in ls
            qs = picksof(t, l)
            scatter!(ax, x[[q.plane for q in qs]], [q.shift for q in qs];
                     color=:transparent, strokecolor=multipletcolour(m), strokewidth=1,
                     markersize=6)
        end
    end
    planes = [i for i in eachindex(s) if !isnothing(linecentre(t, ls, i))]
    scatter!(ax, x[planes], [linecentre(t, ls, i) for i in planes];
             color=multipletcolour(m), markersize=8)
    return fig
end
