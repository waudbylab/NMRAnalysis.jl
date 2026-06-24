function preparestate(expt::Experiment)
    @debug "Preparing state"
    state = Dict{Symbol,Observable}()

    state[:mode] = Observable(:normal) # options = :normal, :renaming, :renamingstart, :moving, :fitting

    # Background-fit control. :fit_generation is bumped whenever the fit inputs change (peaks,
    # radii, fitting toggled off, window closed); an in-flight fit checks it and aborts if it no
    # longer matches the generation it started with. :fit_task holds the current fit task.
    state[:fit_generation] = Observable(0)
    state[:fit_task] = Observable{Union{Task,Nothing}}(nothing)

    state[:total_peaks] = Observable(length(expt.peaks[]))

    state[:current_slice] = Observable(1)

    state[:current_slice_label] = lift(idx -> slicelabel(expt, idx), state[:current_slice])

    state[:current_mask_x] = Observable(expt.specdata.x[1])
    state[:current_mask_y] = Observable(expt.specdata.y[1])
    state[:current_mask_z] = Observable(expt.specdata.mask[][1])
    onany(expt.specdata.mask, state[:current_slice]) do m, idx
        # Notify x/y as well as z: planes may have different axis sizes (e.g. rdc2d's separate
        # spectra), so the contour/heatmap must resolve all three together or it sees a stale
        # axis against a new-size matrix ("Incompatible input axes").
        state[:current_mask_x][] = expt.specdata.x[idx]
        state[:current_mask_y][] = expt.specdata.y[idx]
        state[:current_mask_z][] = m[idx]
    end

    state[:current_spec_x] = Observable(expt.specdata.x[1])
    state[:current_spec_y] = Observable(expt.specdata.y[1])
    state[:current_spec_z] = Observable(expt.specdata.z[1])
    on(state[:current_slice]) do idx
        state[:current_spec_x][] = expt.specdata.x[idx]
        state[:current_spec_y][] = expt.specdata.y[idx]
        state[:current_spec_z][] = expt.specdata.z[idx]
    end

    state[:current_fit_x] = Observable(expt.specdata.x[1])
    state[:current_fit_y] = Observable(expt.specdata.y[1])
    state[:current_fit_z] = Observable(expt.specdata.zfit[][1])
    onany(expt.specdata.zfit, state[:current_slice]) do zfit, idx
        state[:current_fit_x][] = expt.specdata.x[idx]
        state[:current_fit_y][] = expt.specdata.y[idx]
        state[:current_fit_z][] = zfit[idx]
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
        initialpositions = [initialposition(peak)[][idx] for peak in expt]
        positions = [position(peak)[][idx] for peak in expt]
        labels = [peak.label[] for peak in expt]
        touched = [peak.touched[] for peak in expt]
        d = Dict(:initialpositions => initialpositions,
                 :positions => positions,
                 :labels => labels,
                 :touched => touched)
        @debug "current_peaks lift" d maxlog = 10
        d
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
        cols = map(t -> t ? :red : :blue, d[:touched])
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
        notify(state[:labels])
    end

    on(state[:current_peak_idx]) do idx
        @debug "current peak index changed to $idx - updating colours"
        d = state[:current_peaks][]
        cols = map(t -> t ? :red : :blue, d[:touched])
        sizes = fill(15.0, length(cols))
        if idx > 0
            cols[idx] = :lime
            sizes[idx] = 24.0
        end
        state[:peakcolours].val = cols
        state[:initialpeaksizes].val = sizes
        notify(state[:initialpeaksizes])
        notify(state[:peakcolours])
    end

    state[:current_peak_info] = lift(idx -> peakinfotext(expt, idx),
                                     state[:current_peak_idx])

    completestate!(state, expt)

    return state
end

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