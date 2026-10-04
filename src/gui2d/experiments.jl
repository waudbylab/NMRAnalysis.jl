"""
    Experiment

An abstract type representing a type of experimental analysis.

# Implementation Requirements

Concrete subtypes must implement:
- `hasfixedpositions(expt)`: Check if the experiment has fixed peak positions between spectra
- `addpeak!(expt, position)`: Add a peak to the experiment

Expected fields:
- `peaks`: A list of peaks in the experiment
- `specdata`: A SpecData object representing the observed and simulated data plus mask
- `clusters`: An Observable list of clusters of peaks
- `touched`: An Observable list of touched clusters
- `isfitting`: An Observable boolean indicating if real-time fitting is active
- `skipplanes`: Planes excluded from fitting

# Functions handled by the abstract type

- `nslices(expt)`: Get the number of slices in the experiment
- `npeaks(expt)`: Get the number of peaks in the experiment
- `mask!(expt)`, `mask(peaks, expt)`: Calculate peak masks
- `simulate!([z], [peaks], expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the peaks in the experiment (see `fitting.jl`)

"""

# implementations
include("expt-intensitybased.jl")
include("expt-moving.jl")
include("expt-hetnoe.jl")
include("expt-cest.jl")
include("expt-cpmg.jl")
include("expt-ccr.jl")
include("expt-methylccr.jl")

# generic functions

"""
    nslices(expt::Experiment)

Number of spectra in the experiment.
"""
nslices(expt::Experiment) = length(expt.specdata.z)

"""
    npeaks(expt::Experiment) 

Number of peaks currently in the experiment.
"""
npeaks(expt::Experiment) = length(expt.peaks[])

hasfixedpositions(expt::FixedPeakExperiment) = true
hasfixedpositions(expt::MovingPeakExperiment) = false

"""
    primaryparam(expt) -> Symbol

The experiment's primary derived (post-fit) result. It is written first among
the derived columns in `results.csv` and is the default parameter plotted by
`summaryplot`. Defaults to `:amp` (peak amplitude) for experiments that derive
no post-fit parameters (e.g. plain `fit2d`).
"""
primaryparam(::Experiment) = :amp

function setupexptobservables!(expt)
    xres = abs(expt.specdata.x[1][2] - expt.specdata.x[1][1])
    yres = abs(expt.specdata.y[1][2] - expt.specdata.y[1][1])
    expt.xradius[] = clamp(4 * xres, 0.04, 0.1)
    expt.yradius[] = clamp(4 * yres, 0.1, 0.8)
    on(expt.peaks) do _
        @debug "Peaks changed"
        expt.state[][:deferupdates][] && return
        # N.B. do NOT bump :fit_generation here - this observer also fires on the
        # fit's own completion notify(expt.peaks). Genuine user changes reach
        # fit!(expt) (via the touched chain), which bumps the generation itself and
        # supersedes any in-flight fit.
        mask!(expt)
        cluster!(expt)
        return simulate!(expt)
    end
    on(expt.clusters) do _
        @debug "Clusters changed"
        return checktouched!(expt)
    end
    on(expt.touched) do _
        @debug "Cluster touched"
        # touch peaks inside touched clusters
        for i in 1:length(expt.touched[])
            if expt.touched[][i]
                for j in expt.clusters[][i]
                    expt.peaks[][j].touched[] = true
                end
            end
        end
        if expt.isfitting[]
            fit!(expt)
        end
    end
    on(expt.isfitting) do _
        @debug "Fitting changed"
        if expt.isfitting[]
            fit!(expt)
        else
            # fitting toggled off - cancel any in-flight fit (this bypasses fit!,
            # so bump the generation and clear the fitting status ourselves)
            expt.state[][:fit_generation][] += 1
            expt.state[][:mode][] = :normal
        end
    end
    # A change of default radius reaches every peak but those given radii of their own
    # (notify reaches fit!, which supersedes any in-flight fit)
    on(expt.xradius) do r
        for peak in expt.peaks[]
            peak.customradius[] && continue
            peak.xradius.val = r
            peak.touched.val = true
        end
        return notify(expt.peaks)
    end
    on(expt.yradius) do r
        for peak in expt.peaks[]
            peak.customradius[] && continue
            peak.yradius.val = r
            peak.touched.val = true
        end
        return notify(expt.peaks)
    end
    return expt
end

"Mark peak `idx` and the rest of its cluster as needing a refit."
function touchcluster!(expt, idx)
    cluster = findfirst(c -> idx in c, expt.clusters[])
    members = isnothing(cluster) ? [idx] : expt.clusters[][cluster]
    for i in members
        expt.peaks[][i].touched.val = true
    end
    return expt
end

"""
    movepeak!(expt, idx, newpos)

Move peak `idx` to `newpos`. Updates clusters automatically.
"""
function movepeak!(expt, idx, newpos)
    touchcluster!(expt, idx)
    slice = expt.state[][:current_slice][]
    peak = expt.peaks[][idx]
    peak.parameters[:x].initialvalue[][slice] = newpos[1]
    peak.parameters[:y].initialvalue[][slice] = newpos[2]
    # Resample the amplitude at the new position (moving-peak experiments); no-op otherwise.
    reinitialise_amplitude!(expt, peak, slice)
    peak.touched.val = true
    return notify(expt.peaks)
end

# Hook: re-estimate a peak's amplitude in plane `slice` after its position moves. Specialised
# for moving-peak experiments (see expt-moving.jl); a no-op for fixed-position experiments.
reinitialise_amplitude!(::Experiment, peak, slice) = nothing

"""
    deletepeak!(expt, idx)

Delete peak `idx`. Updates clusters automatically.
"""
function deletepeak!(expt, idx)
    touchcluster!(expt, idx)
    deleteat!(expt.peaks[], idx)
    return notify(expt.peaks)
end

# Steps for the Shift+arrow keys, matching the radius sliders
const XRADIUS_STEP = 0.005
const YRADIUS_STEP = 0.02

"""
    bumpradius!(expt, idx, key)

Widen or narrow peak `idx`'s radii by one step, x with the left and right arrows and y with
up and down, giving it radii of its own.
"""
function bumpradius!(expt, idx, key)
    peak = expt.peaks[][idx]
    dx = key == Keyboard.right ? XRADIUS_STEP : key == Keyboard.left ? -XRADIUS_STEP : 0.0
    dy = key == Keyboard.up ? YRADIUS_STEP : key == Keyboard.down ? -YRADIUS_STEP : 0.0
    touchcluster!(expt, idx)
    setradius!(peak, max(peak.xradius[] + dx, XRADIUS_STEP),
               max(peak.yradius[] + dy, YRADIUS_STEP))
    return notify(expt.peaks)
end

"Return peak `idx` to the experiment's default radii."
function resetradius!(expt, idx)
    peak = expt.peaks[][idx]
    peak.customradius[] || return nothing
    touchcluster!(expt, idx)
    setradius!(peak, expt.xradius[], expt.yradius[])
    peak.customradius.val = false
    return notify(expt.peaks)
end

"""
    deleteallpeaks!(expt)

Remove all peaks from the experiment.
"""
function deleteallpeaks!(expt)
    expt.state[][:current_peak_idx][] = 0
    empty!(expt.peaks[])
    return notify(expt.peaks)
end

"""
    batchupdate(f, expt)

Call `f()` with the masking, clustering and fitting that follow each change to the peaks
held back, then do them once. Adding a peak list one peak at a time would otherwise redo
all three, and start and abandon a fit, for every peak.
"""
function batchupdate(f, expt::Experiment)
    deferring = expt.state[][:deferupdates]
    deferring[] = true
    try
        f()
    finally
        deferring[] = false
    end
    return notify(expt.peaks)
end

"""
    checktouched!(expt)

Update which clusters have been modified.
"""
function checktouched!(expt)
    touched = map(expt.clusters[]) do cluster
        return any(j -> expt.peaks[][j].touched[], cluster)
    end
    return expt.touched[] = touched
end

"""
    skipset(expt) -> Set{Int}

Planes excluded from fitting the lineshapes and from the model fitted after them. They are
still displayed, and their amplitudes still measured.
"""
skipset(expt::Experiment) = Set{Int}(expt.skipplanes)

"""
    fit!(expt; warm=false)

Fit every touched cluster in the background, then the post-fits, then notify the peaks once.

A newer fit, a change to the peaks, switching fitting off or closing the window supersedes
this one: it stops at its next step and commits nothing. A cluster that runs past its time
budget keeps its previous values and is marked `:timeout`; `warm=true` continues such
clusters from their last values with a longer budget (see [`continuefit!`](@ref)).
"""
function fit!(expt::Experiment; warm=false)
    state = expt.state[]
    peaks = copy(expt.peaks[])
    todo = [cluster for (cluster, t) in zip(expt.clusters[], expt.touched[]) if t]
    isempty(todo) && return nothing

    mygen = (state[:fit_generation][] += 1)
    budget = warm ? CONTINUE_TIME_BUDGET : FIT_TIME_BUDGET
    clusters = [peaks[cluster] for cluster in todo]

    runfit = function ()
        current() = state[:fit_generation][] == mygen
        try
            state[:mode][] = :fitting
            state[:fitprogress][] = (0, length(clusters))
            statuses = fitclusters!(expt, clusters, mygen, budget; warm)
            (current() && :cancelled ∉ statuses) || return nothing
            for (cluster, status) in zip(clusters, statuses)
                for peak in cluster
                    peak.touched.val = false
                    peak.fitstatus.val = status
                    status == :timeout || postfit!(peak, expt)
                end
            end
            postfitglobal!(expt)
            current() && notify(expt.peaks)
        catch e
            # a background task's error is otherwise never seen
            @error "Fitting failed" exception = (e, catch_backtrace())
        finally
            current() && (state[:mode][] = :normal)
        end
        return nothing
    end

    # A task on the main thread, since it updates observables and so plots; the clusters
    # themselves are fitted on worker threads where there are any (see fitclusters!)
    return state[:fit_task][] = @async runfit()
end

"""
    continuefit!(expt)

Refit the peaks whose last fit was stopped before it finished (see
[`iscontinuable`](@ref)), starting from where it stopped and with a longer time budget.
"""
function continuefit!(expt::Experiment)
    unfinished = filter(iscontinuable, expt.peaks[])
    isempty(unfinished) && return nothing
    foreach(peak -> peak.touched.val = true, unfinished)
    expt.touched.val = map(c -> any(j -> expt.peaks[][j].touched[], c), expt.clusters[])
    return fit!(expt; warm=true)
end

# Additional fitting of a peak following the spectrum fit - by default, none.
postfit!(peak::Peak, expt::Experiment) = (peak.postfitted[] = true)

# Global fitting across every peak following the spectrum fit - by default, none.
postfitglobal!(expt::Experiment) = nothing

"""
    simulate!(expt::Experiment)

Simulate all peaks in the experiment and update specdata.
"""
function simulate!(expt::Experiment)
    z = expt.specdata.zfit.val
    foreach(zi -> fill!(zi, 0), z)
    for peak in expt.peaks[]
        simulate!(z, peak, expt)
    end
    return notify(expt.specdata.zfit)
end

"""
    mask!(expt::Experiment)

Calculate masks for all peaks and update specdata.
"""
function mask!(expt::Experiment)
    z = expt.specdata.mask.val
    for i in eachindex(z)
        if i > 1 && hasfixedpositions(expt) && samegrid(expt, i, i - 1)
            z[i] .= z[i - 1]
        else
            fill!(z[i], false)
            foreach(peak -> maskplane!(z[i], peak, expt, i), expt.peaks[])
        end
    end
    return notify(expt.specdata.mask)
end

"""
    mask(peaks, expt) -> Vector{BitMatrix}

The mask of a cluster of peaks, one matrix per plane. Planes sharing a grid in a
fixed-peak experiment share one matrix.
"""
function mask(peaks::AbstractVector{Peak}, expt::Experiment)
    z = Vector{BitMatrix}(undef, nslices(expt))
    for i in eachindex(z)
        if i > 1 && hasfixedpositions(expt) && samegrid(expt, i, i - 1)
            z[i] = z[i - 1]
        else
            z[i] = falses(size(expt.specdata.z[i]))
            foreach(peak -> maskplane!(z[i], peak, expt, i), peaks)
        end
    end
    return z
end

# The fitting region of `peak` in plane `i`: an ellipse of its radii about where it was
# placed. Specialise this for a different shape of region.
function maskplane!(m, peak::Peak, expt::Experiment, i)
    return maskellipse!(m, expt.specdata.x[i], expt.specdata.y[i],
                        initialposition(peak)[i]..., peak.xradius[], peak.yradius[])
end

function Base.show(io::IO, expt::Experiment)
    return print(io, "$(typeof(expt))($(npeaks(expt)) peaks, $(nslices(expt)) slices)")
end

function Base.show(io::IO, mime::MIME"text/plain", expt::Experiment)
    println(io, "$(typeof(expt))")
    println(io, "  $(npeaks(expt)) peaks")
    println(io, "  $(nslices(expt)) slices")
    println(io, "  $(length(expt.clusters[])) clusters")
    unfinished = count(isunfinished, expt.peaks[])
    unfinished > 0 && println(io, "  $unfinished unfinished fits")
    return print(io,
                 "Use results(expt) for one row per peak, planeresults(expt) for one " *
                 "per plane")
end
