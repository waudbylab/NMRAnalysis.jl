# Reactive state for the interactive GUI. Observables live only here and in `gui.jl`; the
# analysis itself stays pure (`analyse(expt, dataset, regions)`), recomputed live as the
# user moves regions. Cheap derivations use `lift`; nothing expensive is button-gated here
# because integration + a few small fits are fast.

"""
    prepare_state(expt) -> Dict{Symbol,Any}

Build the Observable graph for `expt`: region centres/widths, noise position, current
plane, and the live analysis result with its derived plot data.
"""
function prepare_state(expt::Experiment1D)
    state = Dict{Symbol,Any}()
    state[:expt] = expt
    ds = dataset(expt)
    planes = ds.planes
    state[:planes] = planes
    state[:nplanes] = nplanes(planes)

    # regions as (label, centre, width); width defaults to 0.1 ppm for zero-width regions
    regs0 = map(regions(expt)) do r
        w = r.hi - r.lo
        return (; label=r.label, c=(r.lo + r.hi) / 2, w=(w == 0 ? 0.1 : w))
    end
    state[:regions] = Observable(collect(regs0))
    state[:active] = Observable(1)

    n = ds.noise
    nw = n.hi - n.lo
    state[:noisec] = Observable((n.lo + n.hi) / 2)
    state[:noisew] = Observable(nw == 0 ? 0.2 : nw)

    state[:currentplane] = Observable(1)
    state[:isfitting] = Observable(true)
    state[:outputdir] = Observable("out")

    # current Region objects + dataset (noise applied)
    state[:regionobjs] = lift(state[:regions]) do rs
        return [Region(r.label, r.c - r.w / 2, r.c + r.w / 2) for r in rs]
    end
    state[:noiseobj] = lift(state[:noisec], state[:noisew]) do c, w
        return Region("noise", c - w / 2, c + w / 2)
    end
    state[:dataset] = lift(nr -> Dataset1D(planes, nr), state[:noiseobj])

    # live analysis (independent of the fitting toggle; the toggle only hides fit curves)
    state[:result] = lift(state[:dataset], state[:regionobjs]) do ds_, regs_
        return analyse(expt, ds_, regs_)
    end

    # spectra overlay (static) and current plane
    state[:spectra] = [Point2f.(t.δ, t.y) for t in planes.traces]
    state[:overlay] = _overlay_points(state[:spectra])
    state[:currentspectrum] = lift(i -> state[:spectra][clamp(i, 1, length(state[:spectra]))],
                                   state[:currentplane])

    # active region label and derived plot data
    state[:activelabel] = lift(state[:regions], state[:active]) do rs, i
        return rs[clamp(i, 1, length(rs))].label
    end
    state[:plotdata] = lift(state[:result], state[:activelabel]) do res, lbl
        return result_plotdata(expt, res, lbl)
    end
    state[:plotpoints] = lift(d -> d[1], state[:plotdata])
    state[:ploterrors] = lift(d -> d[2], state[:plotdata])
    state[:plotfit] = lift(state[:plotdata], state[:isfitting]) do d, fitting
        return fitting ? d[3] : Point2f[]
    end
    state[:summary] = lift(res -> summary_text(expt, res), state[:result])

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

# region mutation helpers (copy-and-replace so the Observable fires)
function set_region_center!(state, i, c)
    rs = copy(state[:regions][])
    r = rs[i]
    rs[i] = (; label=r.label, c=c, w=r.w)
    return state[:regions][] = rs
end

function set_region_width!(state, i, w)
    rs = copy(state[:regions][])
    r = rs[i]
    rs[i] = (; label=r.label, c=r.c, w=w)
    return state[:regions][] = rs
end
