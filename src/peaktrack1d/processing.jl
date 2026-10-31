# Time-domain data, and how they become the spectra shown.

"""
    Processing(fids; lb=1.0, zerofill=2, solvent=:auto, phase=:auto)

Time-domain data and how they become the spectra shown: one FID per experiment (1D or
pseudo-2D), exponential line broadening `lb` (Hz, adjustable afterwards), zero-filling,
solvent suppression as in `autoprocess`, and one `(p0, p1)` phase correction per
experiment (degrees, TopSpin's convention) with where it came from:

- `:auto`: `autoprocess`'s automatic phasing
- `:refit`: automatic phasing whose first-order phase was an outlier from the series, so
  its zero-order phase was fitted again with the first-order phase fixed (see
  [`fixfirstorder!`](@ref))
- `:pdata`: PHC0 and PHC1 from the experiment's processed data
- `:given`: passed in, or read from a `processing.csv`

`phase` is `:auto`, `:pdata`, a vector of `(p0, p1)` pairs, or the path of a
`processing.csv` written by an earlier analysis. Phases are found once; changing the line
broadening reprocesses with them, since an exponential window does not change the phase.
"""
struct Processing
    fids::Vector{NMRData}
    phases::Vector{NTuple{2,Float64}}
    phasing::Vector{Symbol}
    zerofill::Int
    solvent::Any
    lb::Observable{Float64}
end

function Processing(fids::AbstractVector{<:NMRData}; lb=1.0, zerofill=2, solvent=:auto,
                    phase=:auto)
    fids = collect(NMRData, fids)
    phases, phasing = initialphases(fids, phase; window=ExponentialWindow(lb), zerofill,
                                    solvent)
    return Processing(fids, phases, phasing, zerofill, solvent, Observable(Float64(lb)))
end

function initialphases(fids, phase::Symbol; kw...)
    phase === :auto && return autophases(fids; kw...)
    phase === :pdata && return pdataphases(fids; kw...)
    throw(ArgumentError("phase must be :auto, :pdata, a vector of (p0, p1) pairs or the " *
                        "path of a processing.csv file, not :$phase"))
end

function initialphases(fids, phases::AbstractVector; kw...)
    length(phases) == length(fids) ||
        throw(DimensionMismatch("$(length(phases)) phases for $(length(fids)) experiments"))
    return [(Float64(p[1]), Float64(p[2])) for p in phases], fill(:given, length(fids))
end

initialphases(fids, path::AbstractString; kw...) = initialphases(fids, readphases(path); kw...)

function autophases(fids; kw...)
    phases = map(enumerate(fids)) do (k, fid)
        @info "Autophasing experiment $k of $(length(fids))"
        return autophase1(fid; kw...)
    end
    phasing = fill(:auto, length(fids))
    fixfirstorder!(phases, phasing, fids; kw...)
    return phases, phasing
end

"The phases `autoprocess` settles on for one FID."
function autophase1(fid; window, zerofill, solvent)
    spec = autoprocess(fid; window, zerofill, solvent, phase=:auto)
    ph = something(metadata(dims(spec, 1), :phase), (p0=0.0, p1=0.0))
    return (Float64(ph.p0), Float64(ph.p1))
end

"""
PHC0 and PHC1 from each experiment's `pdata/1/procs`, falling back to automatic phasing
where there is none. They are checked like automatic phases, since a stored phase is no more
reliable than the processing that left it there.
"""
function pdataphases(fids; kw...)
    phases = NTuple{2,Float64}[]
    phasing = Symbol[]
    for (k, fid) in enumerate(fids)
        procs = joinpath(dirname(string(something(fid[:filename], ""))), "pdata", "1", "procs")
        if isfile(procs)
            dic = loadjdx(procs)
            push!(phases, (Float64(dic[:phc0]), Float64(dic[:phc1])))
            push!(phasing, :pdata)
        else
            @warn "No processed data for experiment $k, so it was phased automatically" procs
            push!(phases, autophase1(fid; kw...))
            push!(phasing, :auto)
        end
    end
    fixfirstorder!(phases, phasing, fids; kw...)
    return phases, phasing
end

"Smallest departure (degrees) from the series median that marks a first-order phase as wrong."
const P1TOLERANCE = 15.0

"""
    fixfirstorder!(phases, phasing, fids; window, zerofill, solvent)

Correct first-order phases that are outliers from the series.

Residual first-order phase comes from the acquisition (the digital filter's group delay,
the pre-scan delay, the excitation pulse) and not from the sample, so across a series
recorded with the same parameters it barely varies, while the zero-order phase follows the
receiver, the tuning and the temperature. A phasing that has settled on a wrong first-order
phase, getting one line right and another wrong, therefore stands out from the series
median. Its zero-order phase is fitted again with the first-order phase fixed at the median.

Outliers are found by [`p1outliers`](@ref).
"""
function fixfirstorder!(phases, phasing, fids; window, zerofill, solvent)
    p1 = last.(phases)
    med, outliers = p1outliers(p1)
    for k in outliers
        spec = autoprocess(fids[k]; window, zerofill, solvent, phase=nothing, real=false)
        est = try
            estimatephase(spec, 1; order=0, p1=med)
        catch err
            err isa NMRToolsError || rethrow()
            nothing
        end
        if isnothing(est)
            @warn "The first-order phase of experiment $k, $(round(p1[k]; digits=1))°, is " *
                  "far from the series median of $(round(med; digits=1))°, but its " *
                  "zero-order phase could not be fitted again, so it was left as it was"
            continue
        end
        @info "The first-order phase of experiment $k, $(round(p1[k]; digits=1))°, is far " *
              "from the series median of $(round(med; digits=1))°, so it was set to the " *
              "median and the zero-order phase fitted again"
        phases[k] = (Float64(est.p0), med)
        phasing[k] = :refit
    end
    return phases
end

"""
    p1outliers(p1) -> (median, indices)

The series median of the first-order phases `p1`, and the indices of those departing from it
by more than five robust standard deviations and by more than `P1TOLERANCE`. None for fewer
than three.
"""
function p1outliers(p1)
    med = median(p1)
    length(p1) < 3 && return med, Int[]
    tolerance = max(P1TOLERANCE, 5 * 1.4826 * median(abs.(p1 .- med)))
    return med, findall(x -> abs(x - med) > tolerance, p1)
end

"FID `k` processed with exponential line broadening `lb` and its stored phases."
function spectrum(p::Processing, k, lb=p.lb[])
    return autoprocess(p.fids[k]; window=ExponentialWindow(lb), zerofill=p.zerofill,
                       solvent=p.solvent, phase=p.phases[k])
end

"The spectra from `p`, one or more planes per experiment."
SpectrumSeries(p::Processing) = SpectrumSeries([spectrum(p, k) for k in eachindex(p.fids)])

"Reprocess every experiment at line broadening `lb`, writing the planes into `s` in place."
function process!(s::SpectrumSeries, p::Processing, lb=p.lb[])
    for k in eachindex(p.fids)
        b = planesof(spectrum(p, k, lb))
        cols = planesfrom(s, k)
        s.intensities[:, cols] .= b.y
        s.noise[cols] .= b.noise
    end
    return s
end

source(p::Processing, k) = string(something(p.fids[k][:filename], ""))
