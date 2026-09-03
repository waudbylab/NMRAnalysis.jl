# Reactive state for the interactive GUI. Observables live only here and in `gui.jl`; the
# analysis itself stays pure (`analyse(expt, dataset, regions)`), recomputed live as the
# user moves/adds/renames/deletes regions or the noise marker. Cheap derivations use
# `lift`; nothing expensive is button-gated here because integration + a few small fits
# are fast.

const ACTIVE_REGION_COLOR = (:orange, 0.4)
const INACTIVE_REGION_COLOR = (:steelblue, 0.25)
const NOISE_COLOR = (:orchid, 0.3)

# Explicit RGBAf conversion, so every colour flowing into the same Observable/plot-array
# is the same concrete type as `Makie.wong_colors()` returns colours without alpha.
const PALETTE = [RGBAf(c.r, c.g, c.b, 1.0) for c in Makie.wong_colors()]

"""Colour for the `i`-th result series, cycling through a fixed palette."""
seriescolor(i) = PALETTE[mod1(i, length(PALETTE))]

"""
    prepare_state(expt) -> Dict{Symbol,Any}

Build the Observable graph for `expt`: regions (position/width/label), the noise marker
position, current spectrum, GUI interaction mode, and the live analysis result with its
derived plot data.
"""
function prepare_state(expt::Experiment1D)
    state = Dict{Symbol,Any}()
    state[:expt] = expt
    ds = dataset(expt)
    planes = ds.planes
    state[:planes] = planes
    state[:nplanes] = nplanes(planes)

    # regions: (label, centre, width) triples. A region with w == 0 is a height.
    regs0 = map(regions(expt)) do r
        w = width(r)
        return (; label=r.label, c=(r.lo + r.hi) / 2, w=(w == 0 ? 0.05 : w))
    end
    state[:regions] = Observable(collect(regs0))
    state[:active] = Observable(isempty(regs0) ? 0 : 1)

    # noise marker: a single position, no independent width (matches whichever region is
    # being reduced - see reduce_region).
    state[:noisec] = Observable(ds.noisecenter)

    state[:currentspectrum_idx] = Observable(1)
    state[:isfitting] = Observable(true)
    state[:outputdir] = Observable("out")

    # interaction mode: :normal, :addingdrag (press-drag-release with 'A'), :renaming
    state[:mode] = Observable(:normal)
    state[:addstart] = Observable(0.0)
    state[:addcurrent] = Observable(0.0)
    state[:oldlabel] = Observable("")

    # current Region objects + dataset (noise position applied)
    state[:regionobjs] = lift(state[:regions]) do rs
        return [Region(r.label, r.c - r.w / 2, r.c + r.w / 2) for r in rs]
    end
    state[:dataset] = lift(nc -> Dataset1D(planes, nc), state[:noisec])

    # live analysis - the Fitting toggle genuinely disables curve-fitting here (see
    # `analyse`/`series_results`' `isfitting`), not just the plot/text display of it.
    #
    # No special-casing for zero regions: `analyse` returns a `Vector{RegionResult}` for
    # every experiment, fitted or not, empty or not, so the Observable's element type -
    # fixed by its first value - is stable whatever the user does.
    state[:result] = lift(state[:dataset], state[:regionobjs], state[:isfitting]) do ds_,
                                                                                      regs_,
                                                                                      fitting
        return analyse(expt, ds_, regs_; isfitting=fitting)
    end

    # spectra overlay (static) and current spectrum
    state[:spectra] = [Point2f.(t.δ, t.y) for t in planes.traces]
    state[:overlay] = _overlay_points(state[:spectra])
    state[:currentspectrum] = lift(state[:currentspectrum_idx]) do i
        return state[:spectra][clamp(i, 1, length(state[:spectra]))]
    end

    # active region label
    state[:activelabel] = lift(state[:regions], state[:active]) do rs, i
        return (i < 1 || i > length(rs)) ? "" : rs[i].label
    end

    # per-series plot data for the active region (one entry per group, e.g. TRACT's
    # TROSY/anti-TROSY, or STD's saturation frequencies), and flattened/coloured arrays
    # for plotting (Makie plot objects want contiguous arrays, not a variable number of
    # separate plot calls).
    state[:seriesdata] = lift(state[:result], state[:activelabel]) do res, lbl
        return result_plotdata(expt, res, lbl)
    end
    state[:flat_points] = lift(sd -> reduce(vcat, (s.points for s in sd); init=Point2f[]),
                               state[:seriesdata])
    state[:flat_point_colors] = lift(state[:seriesdata]) do sd
        return reduce(vcat,
                      (fill(seriescolor(i), length(s.points)) for (i, s) in enumerate(sd));
                      init=RGBAf[])
    end
    state[:flat_errors] = lift(state[:seriesdata]) do sd
        return reduce(vcat, (s.errors for s in sd);
                      init=Tuple{Float64,Float64,Float64}[])
    end
    state[:flat_error_colors] = state[:flat_point_colors]
    state[:flat_fit] = lift(state[:seriesdata], state[:isfitting]) do sd, fitting
        fitting || return Point2f[]
        pts = Point2f[]
        for s in sd
            isempty(s.fitline) && continue
            append!(pts, s.fitline)
            push!(pts, Point2f(NaN32, NaN32))
        end
        return pts
    end
    state[:flat_fit_colors] = lift(state[:seriesdata], state[:isfitting]) do sd, fitting
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

    # live results text for the results panel, split into the raw per-series fit (e.g.
    # A, R) and the second-stage derived summary (e.g. TRACT's τc), so the GUI can show
    # them as separate blocks. Both go blank when the Fitting toggle is off - fitting
    # itself always runs (nothing here skips the actual curve_fit calls; see `isfitting`
    # in reductions/experiments), but showing numbers while "Fitting" reads as off would
    # be misleading.
    state[:resultsheader] = lift(state[:result], state[:activelabel], state[:isfitting]) do res, lbl,
                                                                                             fitting
        fitting || return ""
        return resultsheader(expt, res, lbl)
    end
    state[:secondaryresult] = lift(state[:result], state[:activelabel],
                                   state[:isfitting]) do res, lbl, fitting
        fitting || return ""
        return secondarytext(expt, res, lbl)
    end

    return state
end

function _overlay_points(spectra)
    pts = Point2f[]
    for s in spectra
        append!(pts, s)
        push!(pts, Point2f(NaN32, NaN32))
    end
    return pts
end

# ---- region mutation (copy-and-replace so the Observable fires) ---------------

function set_region_center!(state, i, c)
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    r = rs[i]
    rs[i] = (; label=r.label, c=c, w=r.w)
    return state[:regions][] = rs
end

"""
    set_region_edge!(state, i, fixed, moving)

Move one edge of region `i` to `moving`, keeping the opposite edge at `fixed` - so both
the centre and width update together. Used for edge-dragging (as opposed to
`set_region_center!`, which drags the whole region and keeps its width fixed).
"""
function set_region_edge!(state, i, fixed, moving)
    (1 ≤ i ≤ length(state[:regions][])) || return
    lo, hi = minmax(fixed, moving)
    w = max(hi - lo, 1e-6)
    rs = copy(state[:regions][])
    r = rs[i]
    rs[i] = (; label=r.label, c=(lo + hi) / 2, w=w)
    return state[:regions][] = rs
end

function setactivelabel!(state, label)
    i = state[:active][]
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    r = rs[i]
    rs[i] = (; label=label, c=r.c, w=r.w)
    return state[:regions][] = rs
end

"""Number the default label uniquely against the current region list."""
function nextlabel(state)
    n = length(state[:regions][]) + 1
    existing = Set(r.label for r in state[:regions][])
    label = "region$n"
    while label in existing
        n += 1
        label = "region$n"
    end
    return label
end

"""
    addregion!(state, lo, hi)

Add a new region spanning `lo..hi` (order-independent) and select it. Callers that want
the user to immediately name it (as the GUI's 'A' shortcut does) follow up by entering
rename mode themselves.
"""
function addregion!(state, lo, hi)
    lo, hi = minmax(lo, hi)
    w = max(hi - lo, 1e-6)
    rs = copy(state[:regions][])
    push!(rs, (; label=nextlabel(state), c=(lo + hi) / 2, w=w))
    state[:regions][] = rs
    return state[:active][] = length(rs)
end

"""
    deleteregion!(state, i)

Remove region `i`, if it exists.
"""
function deleteregion!(state, i)
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    deleteat!(rs, i)
    state[:regions][] = rs
    return state[:active][] = clamp(state[:active][], 0, length(rs))
end
