"""
    Experiment

An abstract type representing a type of experimental analysis.

# Implementation Requirements

Concrete subtypes must implement:
- `hasfixedpositions(expt)`: Check if the experiment has fixed peak positions between spectra
- `addpeak!(expt, position)`: Add a peak to the experiment
- `simulate!(z, peak, expt)`: Simulate a single peak
- `mask!(z, peak, expt)`: Get the mask for a single peak

Expected fields:
- `peaks`: A list of peaks in the experiment
- `specdata`: A SpecData object representing the observed and simulated data plus mask
- `clusters`: An Observable list of clusters of peaks
- `touched`: An Observable list of touched clusters
- `isfitting`: An Observable boolean indicating if real-time fitting is active

# Functions handled by the abstract type

- `nslices(expt)`: Get the number of slices in the experiment
- `npeaks(expt)`: Get the number of peaks in the experiment
- `mask!([z], [peaks], expt)`: Calculate peak masks and update internal specdata
- `simulate!([z], [peaks], expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the peaks in the experiment

"""

# implementations
include("expt-intensitybased.jl")
include("expt-moving.jl")
include("expt-hetnoe.jl")
include("expt-PRE.jl")
include("expt-cest.jl")
include("expt-cpmg.jl")
include("expt-ccr.jl")
include("expt-methylccr.jl")

# generic functions

"""Maximum wall-clock time (seconds) allowed for a single cluster fit before it is aborted.

Interactive fitting re-runs on every peak move, so a hard cap keeps the GUI responsive even
if convergence is slow. The residual function checks elapsed time on every iteration.
The budget is applied per-cluster so that a spectrum with many peaks does not exhaust the
budget on early clusters and cancel later ones."""
const FIT_TIME_BUDGET = 30.0

"""Thrown from within a fit's residual function to abort an in-flight fit (because the inputs
changed, the time budget was exceeded, or the user cancelled). Caught silently in `fit!(expt)`."""
struct FitCancelled <: Exception end

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
    on(expt.xradius) do _
        @debug "X radius changed"
        # update xradius for all peaks (notify reaches fit!, which supersedes any
        # in-flight fit)
        for peak in expt.peaks[]
            peak.xradius.val = expt.xradius[]
            peak.touched.val = true
        end
        return notify(expt.peaks)
    end
    on(expt.yradius) do _
        @debug "Y radius changed"
        # update yradius for all peaks (notify reaches fit!, which supersedes any
        # in-flight fit)
        for peak in expt.peaks[]
            peak.yradius.val = expt.yradius[]
            peak.touched.val = true
        end
        return notify(expt.peaks)
    end
    return expt
end

"""
    movepeak!(expt, idx, newpos)

Move peak `idx` to `newpos`. Updates clusters automatically.
"""
function movepeak!(expt, idx, newpos)
    # touch other peaks in the same cluster before moving and notifying
    clusteridx = findfirst(i -> idx in i, expt.clusters[])
    cluster = expt.clusters[][clusteridx]
    for i in cluster
        # this will touch the moved peak as well as any in the same cluster
        expt.peaks[][i].touched.val = true
    end
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
    # touch other peaks in the same cluster before deleting and notifying
    clusteridx = findfirst(i -> idx in i, expt.clusters[])
    cluster = expt.clusters[][clusteridx]
    for i in cluster
        expt.peaks[][i].touched[] = true
    end
    deleteat!(expt.peaks[], idx)
    return notify(expt.peaks)
end

"""
    deleteallpeaks!(expt)

Remove all peaks from the experiment.
"""
function deleteallpeaks!(expt)
    expt.state[][:current_peak_idx][] = 0
    # delete any existing peaks
    for i in length(expt.peaks[]):-1:1
        deletepeak!(expt, i)
    end
end

"""
    checktouched!(expt)

Update which clusters have been modified.
"""
function checktouched!(expt)
    @debug "Checking touched clusters" #maxlog=10
    touched = map(expt.clusters[]) do cluster
        return any([expt.peaks[][j].touched[] for j in cluster])
    end
    return expt.touched[] = touched
end

"""
    fit!(expt::Experiment)

Fit all touched peaks/clusters in the experiment.
"""
function fit!(expt::Experiment)
    @debug "Fitting experiment" #maxlog=10
    anythingchanged = false
    # first check if anything has been touched and needs fitting
    for i in 1:length(expt.clusters[])
        if expt.touched[][i]
            anythingchanged = true
        end
    end
    anythingchanged || return

    # Bump the fit generation. This supersedes any fit already running: that task's residual
    # function will see the new generation on its next iteration and abort via FitCancelled,
    # and its results will not be committed (the notify below is generation-guarded).
    mygen = (expt.state[][:fit_generation][] += 1)

    # Run the fit off the GUI's critical path. With multiple threads we use a real background
    # thread (Threads.@spawn) so the GUI stays fully responsive; single-threaded we fall back to
    # @async, relying on the yield() inside the residual to give the event loop cycles so the fit
    # can still be cancelled. The fitting code only writes peak parameter arrays in place (without
    # notify), so the only visible update is the generation-guarded notify(expt.peaks) at the end.
    runfit = function ()
        try
            expt.state[][:mode][] = :fitting
            # iterate over touched clusters and fit; t0 is reset per-cluster so each gets
            # its own FIT_TIME_BUDGET rather than sharing one budget across the whole spectrum.
            for i in 1:length(expt.clusters[])
                if expt.touched[][i]
                    fit!(expt.clusters[][i], expt, mygen, time())
                    postfit!(expt.clusters[][i], expt)
                end
            end
            postfitglobal!(expt)
            @debug "Fit finished - notifying peaks" #maxlog=10
            # Only commit results if this fit is still the current one - a superseded fit must
            # not push stale results to the display.
            if expt.state[][:fit_generation][] == mygen
                notify(expt.peaks)
            end
        catch e
            e isa FitCancelled || rethrow()
            @debug "Fit cancelled (generation $mygen superseded)"
        finally
            # Reset the status background only if we are still the current fit; a newer fit will
            # manage its own mode. Guarantees the mode resets even on cancellation/error.
            if expt.state[][:fit_generation][] == mygen
                expt.state[][:mode][] = :normal
            end
        end
    end

    if Threads.nthreads() > 1
        expt.state[][:fit_task][] = Threads.@spawn runfit()
    else
        expt.state[][:fit_task][] = @async runfit()
    end
    return expt.state[][:fit_task][]
end

# placeholder functions for additional fitting following spectrum fit
function postfit!(cluster::Vector{Int}, expt::Experiment)
    @debug "Post-fitting cluster $cluster" #maxlog=10
    for i in cluster
        postfit!(expt.peaks[][i], expt)
    end
end

"""Additional fitting of peak following spectrum fit - defaults to no action"""
function postfit!(peak::Peak, expt::Experiment)
    return peak.postfitted[] = true
end

"""Global fitting of entire experiment following spectrum fit - defaults to no action"""
postfitglobal!(expt::Experiment) = nothing

function fit!(cluster::Vector{Int}, expt::Experiment, mygen=expt.state[][:fit_generation][],
              t0=time())
    @debug "Fitting cluster $cluster" #maxlog=10
    peaks = [expt.peaks[][i] for i in cluster]

    # initial parameters
    p0 = pack(peaks, :initial)
    pmin = pack(peaks, :min)
    pmax = pack(peaks, :max)
    # lmfit (levenberg_marquardt) throws if p0 lies outside the bounds, so clamp first
    p0 = clamp.(p0, pmin, pmax)

    # get mask and bounds
    m = mask(cluster, expt)
    xbounds, ybounds = bounds(m)
    smallmask = [m[i][xbounds[i], ybounds[i]] for i in 1:length(m)]
    bigmask = reduce(vcat, vec.(m))
    smallmask = reduce(vcat, vec.(smallmask))
    zobs = reduce(vcat, vec.(expt.specdata.z))[bigmask]
    @debug "zobs" zobs maxlog = 10

    @debug "pre-allocating zsim" maxlog = 10
    zsim = map(1:nslices(expt)) do i
        return similar(expt.specdata.z[i][xbounds[i], ybounds[i]])
    end
    @debug "zsim" zsim maxlog = 10
    zsimm = similar(zobs)
    # create residual function
    function resid(p)
        @debug "resid (start)" maxlog = 10
        # Abort if this fit has been superseded (peaks/radii changed, fitting toggled off, or
        # window closed) or if it has exceeded its time budget. The residual is evaluated every
        # LM iteration, so it is the natural place to hook cancellation.
        (expt.state[][:fit_generation][] != mygen) && throw(FitCancelled())
        (time() - t0 > FIT_TIME_BUDGET) && throw(FitCancelled())
        # Single-threaded: yield so the GUI event loop gets a chance to run (and register a
        # cancellation) between iterations. Harmless under multithreading.
        Threads.nthreads() == 1 && yield()
        unpack!(copy(p), peaks, :value)
        for i in 1:nslices(expt)
            fill!(zsim[i], 0.0)
        end
        simulate!(zsim, cluster, expt, xbounds, ybounds)
        zsimm .= reduce(vcat, vec.(zsim))[smallmask]
        return zobs - zsimm
    end
    @debug "running fit" maxlog = 10
    # Box-constrained, runtime-bounded fit. Bounds keep peak positions within their radius and
    # R2 within sensible limits. Loose tolerances + maxIter are appropriate for interactive
    # re-fitting; :finite (the default) matches the Float64-routed residual - a ForwardDiff
    # Jacobian would be identically zero here. The wall-clock cap is enforced inside `resid` via
    # FIT_TIME_BUDGET.
    sol = LsqFit.lmfit(resid, p0, Float64[];
                       lower=pmin, upper=pmax,
                       autodiff=:finite,
                       maxIter=200, x_tol=1e-4, g_tol=1e-6)
    pfit = coef(sol)
    perr = stderror(sol)
    @debug "fit complete" pfit maxlog = 10
    unpack!(pfit, peaks, :value)
    unpack!(perr, peaks, :uncertainty)

    # untouch peaks in the cluster
    for peak in peaks
        peak.touched.val = false
    end
    # N.B. the update to peaks[] will be notified in parent function once all fitting is complete
end

"""
    simulate!(expt::Experiment)

Simulate all peaks in the experiment and update specdata.
"""
function simulate!(expt::Experiment)
    @debug "Simulating experiment" #maxlog=10
    z = expt.specdata.zfit.val
    for zi in z
        fill!(zi, 0)
    end
    simulate!(z, expt)
    return notify(expt.specdata.zfit)
end

function simulate!(z, expt::Experiment)
    @debug "Simulating experiment with z" maxlog = 10
    for cluster in expt.clusters[]
        simulate!(z, cluster, expt)
    end
end

function simulate!(z, cluster::Vector{Int}, expt::Experiment, xbounds=nothing,
                   ybounds=nothing)
    @debug "Simulating cluster $cluster" maxlog = 10
    for i in cluster
        simulate!(z, expt.peaks[][i], expt, xbounds, ybounds)
    end
end

"""
    mask!(expt::Experiment)

Calculate masks for all peaks and update specdata.
"""
function mask!(expt::Experiment)
    @debug "Masking experiment" #maxlog=10
    z = expt.specdata.mask.val
    for zi in z
        fill!(zi, false)
    end
    mask!(z, expt)
    return notify(expt.specdata.mask)
end

function mask!(z, expt::Experiment)
    @debug "Masking experiment with z" maxlog = 10
    for peak in expt.peaks[]
        mask!(z, peak, expt)
    end
end

function mask(cluster::Vector{Int}, expt::Experiment)
    @debug "Generating cluster mask" maxlog = 10
    z = [similar(zi) for zi in expt.specdata.mask[]]
    for zi in z
        fill!(zi, false)
    end
    for i in cluster
        mask!(z, expt.peaks[][i], expt)
    end
    return z
end

function bounds(mask)
    @debug "Generating cluster bounds from mask" maxlog = 10
    b = map(mask) do m
        # m is a matrix of booleans
        # get projections ix and iy where any value is true
        ix = vec(any(m; dims=2))
        iy = vec(any(m; dims=1))
        return [ix, iy]
    end
    # reshape into list of ix, and list of iy
    ix = [b[i][1] for i in 1:length(b)]
    iy = [b[i][2] for i in 1:length(b)]
    @debug "bounds" sum.(ix) sum.(iy) maxlog = 10

    return ix, iy
end

# generic masking method - can be specialised if needed
function mask!(z, peak::Peak, expt::Experiment)
    @debug "masking peak $(peak.label)" maxlog = 10
    n = length(z)
    for i in 1:n
        x = expt.specdata.x[i]
        y = expt.specdata.y[i]
        maskellipse!(z[i], x, y,
                     initialposition(peak)[][i][1],
                     initialposition(peak)[][i][2],
                     peak.xradius[], peak.yradius[])
    end
end

function Base.show(io::IO, expt::Experiment)
    return print(io, "$(typeof(expt))($(npeaks(expt)) peaks, $(nslices(expt)) slices)")
end

function Base.show(io::IO, mime::MIME"text/plain", expt::Experiment)
    println(io, "$(typeof(expt))")
    println(io, "  $(npeaks(expt)) peaks")
    println(io, "  $(nslices(expt)) slices")
    println(io, "  $(length(expt.clusters[])) clusters")
    return println(io, "  fitting: $(expt.isfitting[])")
end