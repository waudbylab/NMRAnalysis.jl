# How the stack is laid out on screen, and the geometry of the brush. Observables hold the
# view, but nothing here draws.

"""
    StackView(s; gain=2)

How the series is drawn: plane `i` sits on a baseline at `(i - 1) * spacing`, intensities
divided by `norm` (see [`displayscale`](@ref)) are multiplied by `gain`, and plane `i` is
displaced along the shift axis by `offsets[i]`, which is zero unless the spectra are
aligned to the reference.
"""
struct StackView
    spacing::Float64
    gain::Observable{Float64}
    norm::Float64
    offsets::Observable{Vector{Float64}}
    revision::Observable{Int}   # bumped when the data are reprocessed in place
end

function StackView(s::SpectrumSeries; gain=2)
    return StackView(1.0, Observable(Float64(gain)), displayscale(s),
                     Observable(zeros(length(s))), Observable(0))
end

ypos(Δ, g, norm, i, h) = (i - 1) * Δ + g * h / norm
ypos(v::StackView, i, h) = ypos(v.spacing, v.gain[], v.norm, i, h)
baseline(v::StackView, i) = (i - 1) * v.spacing

"Where a measured shift `δ` in plane `i` appears on screen, and the reverse."
shown(v::StackView, i, δ) = δ + v.offsets[][i]
measured(v::StackView, i, x) = x - v.offsets[][i]

function location(v::StackView, q::Pick)
    return Point2f(shown(v, q.plane, q.shift), ypos(v, q.plane, q.height))
end

"""
The plane whose trace, within `w` ppm of `p` on screen, comes closest to `p` vertically:
the plane whose peak the cursor is on, or whose baseline it is nearest.
"""
function nearestrow(s::SpectrumSeries, v::StackView, p, w)
    best, dmin = nothing, Inf
    for i in eachindex(s)
        idx = window(s, measured(v, i, p[1]), w)
        isempty(idx) && continue
        lo, hi = ypos.(Ref(v), i, extrema(view(intensities(s, i), idx)))
        d = lo ≤ p[2] ≤ hi ? 0.0 : min(abs(p[2] - lo), abs(p[2] - hi))
        d < dmin && ((best, dmin) = (i, d))
    end
    return best
end

# ---------------------------------------------------------------- screen geometry

"Pixels per data unit along x and y for the current view of `ax`."
pixelscale(ax::Axis) = Vec2f(widths(ax.scene.viewport[]) ./ widths(ax.finallimits[]))

pixeldistance(a, b, scale) = sqrt(sum(abs2, (a .- b) .* scale))

"""
    hit(t, v, p, scale, r, line=nothing) -> Union{Tuple{Int,Int},Nothing}

Key of the pick whose marker is nearest `p`, if within `r` pixels: among the picks of
`line`, or of every line when `line` is `nothing`.
"""
function hit(t::Tracks, v::StackView, p, scale, r, line=nothing)
    best, dmin = nothing, Float32(r)
    for (k, q) in t.picks
        isnothing(line) || first(k) == line || continue
        d = pixeldistance(location(v, q), p, scale)
        d ≤ dmin && ((best, dmin) = (k, d))
    end
    return best
end

"Is point `m` within `r` pixels of the stroke `pts`?"
covered(pts, m, scale, r) = any(q -> pixeldistance(m, q, scale) ≤ r, pts)

"""
    densify(pts, step, scale) -> Vector{Point2f}

The polyline `pts` with points added so that neighbours are at most `step` pixels apart,
so that a stroke drawn quickly still paints a continuous band.
"""
function densify(pts, step, scale)
    out = Point2f[]
    for (a, b) in zip(pts, Iterators.drop(pts, 1))
        n = max(1, ceil(Int, pixeldistance(a, b, scale) / step))
        append!(out, (a + (b - a) * (k / n) for k in 0:(n - 1)))
    end
    isempty(pts) || push!(out, last(pts))
    return out
end

"""
    painted(band, r, scale, y0, y1) -> Vector{Tuple{Float64,Float64}}

The intervals of displayed shift covered by a brush of radius `r` pixels drawn through the
(densified) points `band`, within the horizontal strip from `y0` to `y1`, merged and in
increasing order. This is the shaded area on screen.
"""
function painted(band, r, scale, y0, y1)
    spans = Tuple{Float64,Float64}[]
    for q in band
        dy = (q[2] < y0 ? y0 - q[2] : q[2] > y1 ? q[2] - y1 : 0.0) * scale[2]   # pixels to the strip
        dy > r && continue
        dx = sqrt(r^2 - dy^2) / scale[1]
        push!(spans, (q[1] - dx, q[1] + dx))
    end
    sort!(spans)
    merged = Tuple{Float64,Float64}[]
    for (lo, hi) in spans
        if !isempty(merged) && lo ≤ last(merged)[2]
            merged[end] = (first(merged[end]), max(last(merged)[2], hi))
        else
            push!(merged, (lo, hi))
        end
    end
    return merged
end

"""
    track(s, v, pts, r, scale; snr=5) -> Vector{Pick}

The peaks painted by a brush of radius `r` pixels drawn through `pts`. Each plane owns the
strip between its baseline and the next one up; in each plane the brush crosses, the tallest
local maximum whose position lies inside the shaded area within that strip is picked,
ignoring maxima below `snr` times the plane's noise. Narrow the brush to keep a strong
neighbour out.
"""
function track(s::SpectrumSeries, v::StackView, pts, r, scale; snr=5)
    picks = Pick[]
    band = densify(pts, r / 2, scale)
    isempty(band) && return picks
    δ = shifts(s)
    for i in eachindex(s)
        spans = painted(band, r, scale, baseline(v, i), baseline(v, i) + v.spacing)
        isempty(spans) && continue
        o = v.offsets[][i]
        inside = d -> any(((lo, hi),) -> lo ≤ d + o ≤ hi, spans)
        y = intensities(s, i)
        threshold = snr * noiselevel(s, i)
        ks = [k
              for k in 2:(length(y) - 1)
              if y[k] ≥ y[k - 1] && y[k] > y[k + 1] && y[k] ≥ threshold && inside(δ[k])]
        isempty(ks) && continue
        d, h = vertex(δ, y, ks[argmax(y[ks])])
        push!(picks, Pick(i, d, h))
    end
    return picks
end

"""
    overlay(t, current, g, offsets, Δ, norm) -> (others, siblings, chosen, bars)

What is drawn over the stack: markers and the segments joining them, split into the current
line, the other lines of its multiplet, and everything else; and for each multiplet of more
than one line, a bar across its lines in each plane where its centre is defined, a tick at
the centre, and a dashed path joining the centres. Segments join picks of a line in
adjacent planes only, so a plane without a pick shows as a gap.
"""
function overlay(t::Tracks, current, g, offsets, Δ, norm)
    others = (points=Point2f[], colours=Int[], segments=Point2f[], segcolours=Int[])
    siblings = (points=Point2f[], segments=Point2f[])
    chosen = (points=Point2f[], segments=Point2f[])
    bars = (segments=Point2f[], colours=Int[], path=Point2f[], pathcolours=Int[])
    am = multipletof(t, current)
    position(q) = Point2f(q.shift + offsets[q.plane], ypos(Δ, g, norm, q.plane, q.height))
    for ((l, i), q) in t.picks
        m = multipletof(t, l)
        isnothing(m) && continue
        layer = l == current ? chosen : m == am ? siblings : others
        push!(layer.points, position(q))
        layer === others && push!(others.colours, m)
        next = get(t.picks, (l, i + 1), nothing)
        isnothing(next) && continue
        push!(layer.segments, position(q), position(next))
        layer === others && push!(others.segcolours, m, m)
    end
    for m in eachindex(t.multiplets)
        ls = orderedlines(t, m)
        length(ls) > 1 || continue
        previous = nothing
        for i in eachindex(offsets)
            c = linecentre(t, ls, i)
            if isnothing(c)
                previous = nothing
                continue
            end
            qs = [t.picks[(l, i)] for l in ls]
            y = maximum(q -> ypos(Δ, g, norm, i, q.height), qs) + 0.08Δ
            lo, hi = extrema(q.shift for q in qs) .+ offsets[i]
            x = c + offsets[i]
            tip = Point2f(x, y + 0.06Δ)
            push!(bars.segments, Point2f(lo, y), Point2f(hi, y), Point2f(x, y), tip)
            append!(bars.colours, (m, m, m, m))
            if !isnothing(previous)
                push!(bars.path, previous, tip)
                push!(bars.pathcolours, m, m)
            end
            previous = tip
        end
    end
    return others, siblings, chosen, bars
end
