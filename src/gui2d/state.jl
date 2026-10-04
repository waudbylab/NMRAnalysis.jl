function preparestate(expt::Experiment)
    @debug "Preparing state"
    state = Dict{Symbol,Observable}()

    # :normal, :renaming, :renamingstart, :moving, :adding, :line, :fitting or :saving
    state[:mode] = Observable(:normal)

    # Background-fit control. :fit_generation is bumped whenever the fit inputs change (peaks,
    # radii, fitting toggled off, window closed); an in-flight fit checks it and aborts if it no
    # longer matches the generation it started with. :fit_task holds the current fit task.
    state[:fit_generation] = Observable(0)
    state[:fit_task] = Observable{Union{Task,Nothing}}(nothing)
    # clusters fitted of those in the current fit
    state[:fitprogress] = Observable((0, 0))
    # true while a batch of changes is made (see batchupdate)
    state[:deferupdates] = Observable(false)

    state[:total_peaks] = Observable(length(expt.peaks[]))

    # Output folder, relative to the working directory, typed into the box beside Save.
    state[:outputdir] = Observable("out")

    state[:current_slice] = Observable(1)

    state[:current_slice_label] = lift(state[:current_slice]) do idx
        return slicelabel(expt, idx) * (idx in skipset(expt) ? " [skipped]" : "")
    end

    state[:current_mask_x] = Observable(expt.specdata.x[1])
    state[:current_mask_y] = Observable(expt.specdata.y[1])
    # The displayed mask and fit are copies, never the live planes. A plot input given a new
    # array is compared with the last by value (ComputePipeline's is_same), so an input that
    # aliased a live plane, updated in place, would always compare equal and never redraw.
    state[:current_mask_z] = Observable(copy(expt.specdata.mask[][1]))
    onany(expt.specdata.mask, state[:current_slice]) do m, idx
        # Notify x/y as well as z: planes may have different axis sizes (e.g. rdc2d's separate
        # spectra), so the contour/heatmap must resolve all three together or it sees a stale
        # axis against a new-size matrix ("Incompatible input axes").
        state[:current_mask_x][] = expt.specdata.x[idx]
        state[:current_mask_y][] = expt.specdata.y[idx]
        return state[:current_mask_z][] = copy(m[idx])
    end

    state[:current_spec_x] = Observable(expt.specdata.x[1])
    state[:current_spec_y] = Observable(expt.specdata.y[1])
    state[:current_spec_z] = Observable(expt.specdata.z[1])
    on(state[:current_slice]) do idx
        state[:current_spec_x][] = expt.specdata.x[idx]
        state[:current_spec_y][] = expt.specdata.y[idx]
        return state[:current_spec_z][] = expt.specdata.z[idx]
    end

    state[:current_fit_x] = Observable(expt.specdata.x[1])
    state[:current_fit_y] = Observable(expt.specdata.y[1])
    state[:current_fit_z] = Observable(copy(expt.specdata.zfit[][1]))
    onany(expt.specdata.zfit, state[:current_slice]) do zfit, idx
        state[:current_fit_x][] = expt.specdata.x[idx]
        state[:current_fit_y][] = expt.specdata.y[idx]
        return state[:current_fit_z][] = copy(zfit[idx])
    end

    state[:current_peak_idx] = Observable(0)
    state[:current_peak] = Observable{Union{Peak,Nothing}}(nothing)
    map!(state[:current_peak], expt.peaks, state[:current_peak_idx]) do peaks, idx
        if idx > 0 && idx <= length(peaks)
            peaks[idx]
        else
            nothing
        end
    end

    state[:current_peaks] = lift(expt.peaks, state[:current_slice]) do expt, idx
        initialpositions = [initialposition(peak)[idx] for peak in expt]
        positions = [position(peak)[idx] for peak in expt]
        labels = [peak.label[] for peak in expt]
        colours = [peakcolour(peak) for peak in expt]
        d = Dict(:initialpositions => initialpositions,
                 :positions => positions,
                 :labels => labels,
                 :colours => colours)
        @debug "current_peaks lift" d maxlog = 10
        return d
    end
    state[:initialpositions] = Observable{Vector{Point2f}}([])
    state[:positions] = Observable{Vector{Point2f}}([])
    state[:labels] = Observable{Vector{String}}([])
    state[:oldlabel] = Observable("")
    state[:peakcolours] = Observable{Vector{Symbol}}([])
    # Per-peak handle sizes; the selected/moused-over peak is enlarged for grab feedback. Kept
    # in lockstep with peakcolours (set .val before positions notify) so lengths always match.
    state[:initialpeaksizes] = Observable{Vector{Float64}}([])
    on(state[:current_peaks]) do d
        # onany(state[:current_peaks], state[:current_peak_idx]) do d, idx
        @debug "current peaks changed"
        idx = state[:current_peak_idx][]
        cols = copy(d[:colours])
        sizes = fill(15.0, length(cols))
        if idx > 0
            cols[idx] = :lime
            sizes[idx] = 24.0
        end
        state[:peakcolours].val = cols
        state[:initialpeaksizes].val = sizes
        state[:labels].val = d[:labels]
        state[:positions][] = d[:positions]
        state[:initialpositions][] = d[:initialpositions]
        notify(state[:peakcolours])
        notify(state[:initialpeaksizes])
        return notify(state[:labels])
    end

    on(state[:current_peak_idx]) do idx
        @debug "current peak index changed to $idx - updating colours"
        d = state[:current_peaks][]
        cols = copy(d[:colours])
        sizes = fill(15.0, length(cols))
        if idx > 0
            cols[idx] = :lime
            sizes[idx] = 24.0
        end
        state[:peakcolours].val = cols
        state[:initialpeaksizes].val = sizes
        notify(state[:initialpeaksizes])
        return notify(state[:peakcolours])
    end

    # Depends on expt.peaks as well as the selection, so derived results (e.g. a titration's
    # global Kd) refresh in the panel when a fit completes, not only on reselection.
    state[:current_peak_info] = lift(state[:current_peak_idx], state[:mode],
                                     expt.peaks) do idx, mode, _peaks
        if mode == :adding
            "Adding peak\n\n(a) mark this plane\n(space) fill remaining planes\n(esc) cancel"
        elseif mode == :line
            "Marking a line\n\nMove to its other end, then release (L) if held,\n" *
            "or press (L) again\n(esc) cancel"
        else
            peakinfotext(expt, idx) * radiusnote(expt, idx) * statusnote(expt, idx)
        end
    end

    completestate!(state, expt)

    return state
end

"""
    peakcolour(peak) -> Symbol

Red while a peak awaits fitting, orange where its last fit stopped short of a converged
optimum, blue once fitted.
"""
function peakcolour(peak::Peak)
    peak.touched[] && return :red
    isunfinished(peak) && return :darkorange
    return :blue
end

const STATUS_NOTES = Dict(:maxiter => "it stopped at the iteration limit",
                          :timeout =>
                              "it stopped at the time limit, so the values shown " *
                              "are from before it",
                          :bound =>
                              "a position or linewidth is at its limit. Move the " *
                              "peak, or widen its radius with Shift+arrows")

"The radii of peak `idx` for the info panel, where it has radii of its own."
function radiusnote(expt, idx)
    (idx > 0 && expt.peaks[][idx].customradius[]) || return ""
    peak = expt.peaks[][idx]
    rx = round(peak.xradius[]; digits=3)
    ry = round(peak.yradius[]; digits=2)
    return "\n\nOwn radii: $rx × $ry ppm (press = to reset)"
end

"A warning for the info panel where peak `idx`'s last fit did not converge."
function statusnote(expt, idx)
    (idx > 0 && isunfinished(expt.peaks[][idx])) || return ""
    status = fitstatus(expt.peaks[][idx])
    note = "\n\nFit unfinished: $(STATUS_NOTES[status])."
    return iscontinuable(expt.peaks[][idx]) ?
           note * "\nPress (C) to continue fitting, or Shift+C until it converges." :
           note
end

"""
    initialslice(expt) -> Int

The plane the window opens on: the first, unless the experiment says otherwise.
"""
initialslice(::Experiment) = 1

"Generic label for spectum slices"
slicelabel(expt::Experiment, idx) = "Slice $idx of $(nslices(expt))"

"Generic handler for peak info text"
function peakinfotext(expt::Experiment, idx)
    if idx > 0
        "Peak $idx"
    else
        "No peak selected"
    end
end
