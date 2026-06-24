# Adapters converting NMRData into the pure `Trace`/`Planes`/`Dataset1D` types.
# These are the only part of the analysis layer that touches NMRData; everything
# downstream operates on plain vectors and is testable headless.

"""
    traces_from_spec(spec) -> Vector{Trace}

Extract one `Trace` per plane from a 2D NMRData (dimension 1 = chemical shift,
dimension 2 = planes).
"""
function traces_from_spec(spec)
    δ = collect(data(spec, F1Dim))
    Y = data(spec)
    return [Trace(δ, collect(Y[:, i])) for i in 1:size(Y, 2)]
end

"""
    default_noise_region(spec; frac=0.9, widthfrac=0.05) -> Region

A fallback noise region: a window of width `widthfrac` of the spectral span, centred at
`frac` of the way across the first axis. The GUI lets the user reposition it.
"""
function default_noise_region(spec; frac=0.9, widthfrac=0.05)
    δ = collect(data(spec, F1Dim))
    lo, hi = extrema(δ)
    span = hi - lo
    centre = lo + frac * span
    half = 0.5 * widthfrac * span
    return Region("noise", centre - half, centre + half)
end

"""
    dataset_from_spec(spec, vars; noise=default_noise_region(spec)) -> Dataset1D

Build a `Dataset1D` from a pseudo-2D `spec` and a vector of per-plane variable
`NamedTuple`s (`length(vars) == number of planes`).
"""
function dataset_from_spec(spec, vars::AbstractVector{<:NamedTuple};
                           noise::Region=default_noise_region(spec))
    return Dataset1D(Planes(traces_from_spec(spec), collect(vars)), noise)
end

# ---- relaxation ---------------------------------------------------------------

"""
    load_relaxation(spec; ir=false, tau=nothing, kwargs...) -> RelaxationExperiment

Build a relaxation experiment from a pseudo-2D `spec`. Relaxation delays are taken from
`tau`, or from the `vdlist` in the acquisition parameters when not supplied.
"""
function load_relaxation(spec; ir::Bool=false, tau=nothing, kwargs...)
    times = isnothing(tau) ? acqus(spec, :vdlist) : tau
    vars = [(; time=Float64(t)) for t in times]
    ds = dataset_from_spec(spec, vars)
    return RelaxationExperiment(ds; ir, kwargs...)
end

# ---- TRACT --------------------------------------------------------------------

"""
    load_tract(trosy, antitrosy; tau=nothing, kwargs...) -> TractExperiment

Build a TRACT experiment from TROSY and anti-TROSY pseudo-2D spectra. The two spectra
are concatenated into one dataset tagged by `which ∈ {:trosy, :anti}`. The ¹⁵N Larmor
frequency and cross-correlation prefactor are derived from the TROSY acquisition
parameters.
"""
function load_tract(trosy, antitrosy; tau=nothing, kwargs...)
    ttau = isnothing(tau) ? acqus(trosy, :vdlist) : tau
    atau = isnothing(tau) ? acqus(antitrosy, :vdlist) : tau

    ttraces = traces_from_spec(trosy)
    atraces = traces_from_spec(antitrosy)
    traces = vcat(ttraces, atraces)
    vars = vcat([(; time=Float64(t), which=:trosy) for t in ttau],
                [(; time=Float64(t), which=:anti) for t in atau])

    γH = 2.6752218744e8
    B0 = 2π * acqus(trosy, :bf1) / γH
    ωN = 2π * acqus(trosy, :bf3)
    f = tract_f(; B0)

    ds = Dataset1D(Planes(traces, vars), default_noise_region(trosy))
    return TractExperiment(ds; ωN, f, kwargs...)
end

# ---- nutation calibration -----------------------------------------------------

"""
    load_nutation(spec; durations=nothing, phase=:sine, kwargs...) -> NutationExperiment

Build a nutation calibration experiment from a pseudo-2D `spec` arrayed over pulse
duration. Durations are taken from `durations`, or from the `:calibration`/`:duration`
annotation when not supplied.
"""
function load_nutation(spec; durations=nothing, phase::Symbol=:sine, kwargs...)
    t = isnothing(durations) ? annotations(spec, :calibration, :duration) : durations
    vars = [(; duration=Float64(d)) for d in t]
    ds = dataset_from_spec(spec, vars)
    return NutationExperiment(ds; phase, kwargs...)
end

# ---- STD ----------------------------------------------------------------------

"""
    load_std(spec, sat, tsat; reference=:reference, excess=1.0, regions) -> STDExperiment

Build an STD experiment from a pseudo-2D `spec` whose planes are described by parallel
vectors `sat` (saturation condition per plane) and `tsat` (saturation time per plane).
"""
function load_std(spec, sat::AbstractVector, tsat::AbstractVector; reference=:reference,
             excess::Real=1.0, regions)
    vars = [(; sat=sat[i], tsat=Float64(tsat[i])) for i in eachindex(sat)]
    ds = dataset_from_spec(spec, vars)
    return STDExperiment(ds; regions, reference, excess)
end

# ---- kinetics -----------------------------------------------------------------

"""
    load_kinetics(spec, times; run=nothing, regions, model=NoFitting()) -> KineticsExperiment

Build a kinetics experiment from a pseudo-2D `spec` arrayed over `times`, optionally
tagged by `run` (a per-plane run identifier) for multiple time series.
"""
function load_kinetics(spec, times::AbstractVector; run=nothing, regions,
                  model::SeriesModel=NoFitting())
    vars = if isnothing(run)
        [(; time=Float64(times[i])) for i in eachindex(times)]
    else
        [(; time=Float64(times[i]), run=run[i]) for i in eachindex(times)]
    end
    ds = dataset_from_spec(spec, vars)
    return KineticsExperiment(ds; regions, model)
end
