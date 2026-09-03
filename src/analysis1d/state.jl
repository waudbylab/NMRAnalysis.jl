# Reactive state for the interactive GUI. Observables live only here and in `gui.jl`; the
# analysis itself stays pure (`analyse(expt, dataset, regions)`), recomputed live as the
# user moves/adds/renames/deletes regions or the noise marker. Cheap derivations use
# `lift`; nothing expensive is button-gated here because integration + a few small fits
# are fast.

const ACTIVE_REGION_COLOR = (:orange, 0.4)
const INACTIVE_REGION_COLOR = (:steelblue, 0.25)
const NOISE_COLOR = (:orchid, 0.3)

"Width (ppm) a zero-width region is opened at in the GUI, where it would otherwise be
invisible and impossible to grab."
const MIN_GUI_REGION_WIDTH = 0.05

"""
    preparestate(expt) -> Dict{Symbol,Any}

Build the Observable graph for `expt`: regions (position/width/label), the noise marker
position, current spectrum, GUI interaction mode, and the live analysis result with its
derived plot data.
"""
function preparestate(expt::Experiment1D)
    state = Dict{Symbol,Any}()
    state[:expt] = expt
    ds = dataset(expt)
    planes = ds.planes
    state[:planes] = planes
    state[:nplanes] = nplanes(planes)

    # The GUI holds `Region`s directly - the same type the analysis layer takes - so there
    # is no second representation to keep in step. Mutation is copy-and-replace, which is
    # what makes the Observable fire.
    #
    # A zero-width region (`Region`'s "height" form) is opened at a visible width: it can
    # be neither seen nor grabbed otherwise. Height semantics still apply everywhere the
    # GUI isn't involved.
    state[:regions] = Observable([width(r) == 0 ? setwidth(r, MIN_GUI_REGION_WIDTH) : r
                                  for r in regions(expt)])
    state[:active] = Observable(isempty(state[:regions][]) ? 0 : 1)

    # noise marker: a single position, no independent width (matches whichever region is
    # being reduced - see reduceregion).
    state[:noisec] = Observable(ds.noisecenter)

    state[:currentspectrumidx] = Observable(1)
    state[:isfitting] = Observable(true)
    state[:outputdir] = Observable("out")

    # interaction mode: :normal, :addingdrag (press-drag-release with 'A'), :renaming
    state[:mode] = Observable(:normal)
    state[:addstart] = Observable(0.0)
    state[:addcurrent] = Observable(0.0)
    state[:oldlabel] = Observable("")

    # current Region objects + dataset (noise position applied)
    state[:dataset] = lift(nc -> Dataset1D(planes, nc, ds.label), state[:noisec])

    # live analysis - the Fitting toggle genuinely disables curve-fitting here (see
    # `analyse`/`seriesresults`' `isfitting`), not just the plot/text display of it.
    #
    # No special-casing for zero regions: `analyse` returns a `Vector{RegionResult}` for
    # every experiment, fitted or not, empty or not, so the Observable's element type -
    # fixed by its first value - is stable whatever the user does.
    state[:result] = lift(state[:dataset], state[:regions], state[:isfitting]) do ds_,
                                                                                      regs_,
                                                                                      fitting
        return analyse(expt, ds_, regs_; isfitting=fitting)
    end

    # spectra overlay (static) and current spectrum
    state[:spectra] = [Point2f.(t.δ, t.y) for t in planes.traces]
    state[:overlay] = overlaypoints(state[:spectra])
    state[:currentspectrum] = lift(state[:currentspectrumidx]) do i
        return state[:spectra][clamp(i, 1, length(state[:spectra]))]
    end

    # active region label
    state[:activelabel] = lift(state[:regions], state[:active]) do rs, i
        return (i < 1 || i > length(rs)) ? "" : rs[i].label
    end

    # Result-panel Observables are built by the experiment's visualisation strategy, not
    # here: what the panel needs depends on how it draws (see `ResultVisualisation`).
    completeresultstate!(state, expt)

    # live results text for the results panel, split into the fit's own parameters (e.g.
    # A, R) and the derived quantities post-fitted from them (e.g. TRACT's τc), so the GUI
    # can show them as separate blocks. Both go blank when the Fitting toggle is off -
    # which genuinely skips the `curve_fit` calls (see `isfitting` in `seriesresults`),
    # so there would be nothing to show in any case.
    #
    # Both are rich text (bold group/derived headers over plain parameter lines), not
    # plain String, so the blank branch below returns `rich()` rather than "" - an
    # Observable's element type is fixed by its first value, and a later String wouldn't
    # convert to the RichText the non-blank branch produces (the same hazard `RegionResult`
    # was introduced to avoid).
    state[:resultsheader] = lift(state[:result], state[:activelabel], state[:isfitting]) do res, lbl,
                                                                                             fitting
        fitting || return rich()
        return resultsheader(expt, res, lbl)
    end
    state[:secondaryresult] = lift(state[:result], state[:activelabel],
                                   state[:isfitting]) do res, lbl, fitting
        fitting || return rich()
        return secondarytext(expt, res, lbl)
    end

    return state
end

function overlaypoints(spectra)
    pts = Point2f[]
    for s in spectra
        append!(pts, s)
        push!(pts, Point2f(NaN32, NaN32))
    end
    return pts
end

# ---- region mutation (copy-and-replace so the Observable fires) ---------------

function setregioncentre!(state, i, c)
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    rs[i] = recentre(rs[i], c)
    return state[:regions][] = rs
end

"""
    setregionedge!(state, i, fixed, moving)

Move one edge of region `i` to `moving`, keeping the opposite edge at `fixed` - so both
the centre and width update together. Used for edge-dragging (as opposed to
`setregioncentre!`, which drags the whole region and keeps its width fixed).
"""
function setregionedge!(state, i, fixed, moving)
    (1 ≤ i ≤ length(state[:regions][])) || return
    lo, hi = minmax(fixed, moving)
    rs = copy(state[:regions][])
    rs[i] = Region(rs[i].label, lo, max(hi, lo + 1e-6))
    return state[:regions][] = rs
end

"""
    setregionwidth!(state, i, w)

Resize region `i` to width `w` about its current centre (as opposed to
`setregionedge!`, which changes the centre along with the width). Used by
shift+scroll resizing.
"""
function setregionwidth!(state, i, w)
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    rs[i] = setwidth(rs[i], max(w, 1e-6))
    return state[:regions][] = rs
end

function setactivelabel!(state, label)
    i = state[:active][]
    (1 ≤ i ≤ length(state[:regions][])) || return
    rs = copy(state[:regions][])
    rs[i] = Region(label, rs[i].lo, rs[i].hi)
    return state[:regions][] = rs
end

"""Number the default label uniquely against the current region list. "peak" rather than
"region", so the heading built from it ("Region: peak2") doesn't repeat "region" - and to
match the vocabulary GUI2D uses for the same idea."""
function nextlabel(state)
    n = length(state[:regions][]) + 1
    existing = Set(r.label for r in state[:regions][])
    label = "peak$n"
    while label in existing
        n += 1
        label = "peak$n"
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
    rs = copy(state[:regions][])
    push!(rs, Region(nextlabel(state), lo, max(hi, lo + 1e-6)))
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
