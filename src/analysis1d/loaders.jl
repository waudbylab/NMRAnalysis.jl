# Top-level entry points. Each loads the data, builds the experiment, and opens the
# interactive GUI - matching the behaviour of the routines they replace, and letting the
# analysis-dispatch registry call them directly. Passing an `integration` triple
# `(; peakppm, noiseppm, ppmwidth)` skips the GUI and analyses that region directly, so a
# previously-chosen region can be replayed from a script. These are the only functions in
# the module that touch NMRData; everything downstream works on plain vectors.

_spec(x::AbstractString) = loadnmr(String(x))
_spec(x::Integer) = loadnmr(string(x))
_spec(x) = x

"""Annotation lookup returning `nothing` rather than throwing when absent."""
function _ann(spec, keys...)
    try
        return annotations(spec, keys...)
    catch
        return nothing
    end
end

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
    relaxation1d(spec; ir=nothing, tau=nothing, regions=nothing, integration=nothing)

Analyse a 1D relaxation experiment (T1/T2, or inversion recovery) from a pseudo-2D `spec`
(a path or an `NMRData`).

Relaxation delays come from `tau`, else the `relaxation.duration` annotation, else the
`vdlist`. The model defaults to the `relaxation.model` annotation
(`"inversion_recovery"` selects a recovery fit); pass `ir=true`/`false` to override.
"""
function relaxation1d(spec; ir=nothing, tau=nothing, regions=nothing, integration=nothing)
    spec = _spec(spec)
    times = something(tau, _ann(spec, :relaxation, :duration), acqus(spec, :vdlist))
    isrecovery = isnothing(ir) ?
                 (_ann(spec, :relaxation, :model) == "inversion_recovery") : ir
    ds = dataset_from_spec(spec, [(; time=Float64(t)) for t in times])
    expt = isnothing(regions) ? RelaxationExperiment(ds; ir=isrecovery) :
           RelaxationExperiment(ds; ir=isrecovery, regions)
    return run1d(expt; integration)
end

# ---- TRACT --------------------------------------------------------------------

"""
    tract(trosy, antitrosy; tau=nothing, regions=nothing, integration=nothing)

Analyse a TRACT pair, deriving τc from the TROSY / anti-TROSY relaxation-rate difference.
The two spectra are combined into one dataset tagged by `which ∈ {:trosy, :anti}`, sharing
a single integration region.
"""
function tract(trosy, antitrosy; tau=nothing, regions=nothing, integration=nothing)
    trosy, antitrosy = _spec(trosy), _spec(antitrosy)
    ttau = isnothing(tau) ? acqus(trosy, :vdlist) : tau
    atau = isnothing(tau) ? acqus(antitrosy, :vdlist) : tau

    traces = vcat(traces_from_spec(trosy), traces_from_spec(antitrosy))
    vars = vcat([(; time=Float64(t), which=:trosy) for t in ttau],
                [(; time=Float64(t), which=:anti) for t in atau])

    γH = 2.6752218744e8
    B0 = 2π * acqus(trosy, :bf1) / γH
    ωN = 2π * acqus(trosy, :bf3)
    f = tract_f(; B0)

    ds = Dataset1D(Planes(traces, vars), default_noise_region(trosy))
    expt = isnothing(regions) ? TractExperiment(ds; ωN, f) :
           TractExperiment(ds; ωN, f, regions)
    return run1d(expt; integration)
end

# ---- nutation calibration -----------------------------------------------------

"""
    calibration1d(spec; durations=nothing, phase=nothing, regions=nothing, integration=nothing)

Analyse a 1D nutation calibration, reporting the nutation frequency, 90° pulse length and
B₁ inhomogeneity. Pulse durations come from `durations`, else the `calibration.duration`
annotation. The modulation defaults to the `calibration.model` annotation
(`"cosine_modulated"` selects a cosine); pass `phase=:sine`/`:cosine` to override.
"""
function calibration1d(spec; durations=nothing, phase=nothing, regions=nothing,
                       integration=nothing)
    spec = _spec(spec)
    t = something(durations, _ann(spec, :calibration, :duration))
    ph = isnothing(phase) ?
         (_ann(spec, :calibration, :model) == "cosine_modulated" ? :cosine : :sine) : phase
    ds = dataset_from_spec(spec, [(; duration=Float64(d)) for d in t])
    expt = isnothing(regions) ? NutationExperiment(ds; phase=ph) :
           NutationExperiment(ds; phase=ph, regions)
    return run1d(expt; integration)
end

# ---- diffusion ----------------------------------------------------------------

"""
    diffusion1d(spec, gradients; coherence=SQ(H1), δ=nothing, Δ=nothing, σ=nothing,
                Gmax=0.55, regions=nothing, integration=nothing)

Analyse a diffusion experiment, fitting the Stejskal–Tanner equation to extract the
diffusion coefficient (and the hydrodynamic radius where the solvent and temperature are
known).

`gradients` is the vector of relative gradient strengths (0–1), one per plane. The
gradient pulse length `δ` (`2·p30`), diffusion delay `Δ` (`d20`) and shape factor `σ`
(from `gpnam6`) default to the values in the acquisition parameters.
"""
function diffusion1d(spec, gradients; coherence=SQ(H1), δ=nothing, Δ=nothing, σ=nothing,
                     Gmax=0.55, regions=nothing, integration=nothing)
    spec = _spec(spec)
    γ = gyromagneticratio(coherence)
    δ = isnothing(δ) ? acqus(spec, :p, 30) * 2.0 : δ
    Δ = isnothing(Δ) ? acqus(spec, :d, 20) : Δ
    σ = isnothing(σ) ? _shapefactor(acqus(spec, :gpnam, 6)) : σ

    temp = acqus(spec, :te)
    solvent = _solvent(acqus(spec, :solvent))

    ds = dataset_from_spec(spec, [(; gradient=Float64(g)) for g in gradients])
    kw = (; γ, δ, Δ, σ, Gmax, temp, solvent)
    expt = isnothing(regions) ? DiffusionExperiment(ds; kw...) :
           DiffusionExperiment(ds; kw..., regions)
    return run1d(expt; integration)
end

"""Gradient shape factor from the Bruker shape name (as in the legacy routine)."""
function _shapefactor(gpnam)
    length(gpnam) ≥ 4 || return 1.0
    gpnam[1:4] == "SMSQ" && return 0.9
    gpnam[1:4] == "SINE" && return 0.6366
    return 1.0
end

_solvent(s) = s == "D2O" ? :d2o : (s == "H2O+D2O" ? :h2o : nothing)

# ---- STD ----------------------------------------------------------------------

"""
    std1d(spec, sat, tsat; regions, reference=:reference, excess=1.0, integration=nothing)

Analyse an STD experiment. `sat` and `tsat` are per-plane vectors giving the saturation
condition (one of which is the `reference`) and the saturation time. Reports STD fractions
per region, buildup curves where several saturation times are present, and an epitope map.
"""
function std1d(spec, sat::AbstractVector, tsat::AbstractVector; regions,
               reference=:reference, excess::Real=1.0, integration=nothing)
    spec = _spec(spec)
    vars = [(; sat=sat[i], tsat=Float64(tsat[i])) for i in eachindex(sat)]
    ds = dataset_from_spec(spec, vars)
    return run1d(STDExperiment(ds; regions, reference, excess); integration)
end

# ---- kinetics -----------------------------------------------------------------

"""
    kinetics1d(spec, times; run=nothing, regions, model=NoFitting(), integration=nothing)

Track integrated intensity of one or more named regions against `times`, optionally tagged
by `run` for several time series.
"""
function kinetics1d(spec, times::AbstractVector; run=nothing, regions,
                    model::SeriesModel=NoFitting(), integration=nothing)
    spec = _spec(spec)
    vars = if isnothing(run)
        [(; time=Float64(times[i])) for i in eachindex(times)]
    else
        [(; time=Float64(times[i]), run=run[i]) for i in eachindex(times)]
    end
    ds = dataset_from_spec(spec, vars)
    return run1d(KineticsExperiment(ds; regions, model); integration)
end
