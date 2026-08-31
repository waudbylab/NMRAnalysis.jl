# Methyl cross-correlated relaxation (CCR) analysis.
#
# Two pseudo-3D series — a *buildup* (Iₐ) and a *decay* (I_b) — are recorded as a
# function of a relaxation delay T. For every methyl peak the intensity ratio
# |Iₐ/I_b| is fitted against T to eq 7 of the source paper:
#
#   |Iₐ/I_b| = C·η·tanh(√(η²+δ²)·T) / (√(η²+δ²) − δ·tanh(√(η²+δ²)·T))
#
# with C = 3/4 (triple-quantum) or 1/2 (double-quantum), η the cross-correlated
# relaxation rate (s⁻¹) and δ (< 0) a coupling term. The fitted η is converted to
# the methyl order parameter × tumbling time S²τc (eq 1, ideal methyl geometry).
#
# Implementation: the buildup (N planes) and decay (N planes) datasets are
# concatenated into a single 2N-plane SpecData, so the whole IntensityExperiment
# peak-fitting / GUI / results pipeline is reused unchanged. The CCR-specific
# behaviour lives in `MethylCCRModel`-dispatched methods (`postfit!`,
# `get_model_data`, `setup_post_parameters!`, `model_parameter_text`).

"""
    METHYL_K

Prefactor K (s⁻²) in eq 1 relating η to S²τc for an ideal methyl group, so that
`S²τc = η / K`. Assumes the H–H vector is perpendicular to the methyl 3-fold axis
(θ = 90°, [P₂(cosθ)]² = 1/4) and r_HH = 1.813 Å:

    K = (9/40)·(μ₀/4π)²·γ_H⁴·ℏ² / r_HH⁶  ≈ 3.61×10⁹ s⁻²

giving `S²τc (ns) ≈ 0.277·η` for η in s⁻¹.
"""
const METHYL_K = let
    μ0_4π = 1.0e-7          # μ₀/4π  / T m A⁻¹
    γH = 2.6752218744e8    # ¹H gyromagnetic ratio / rad s⁻¹ T⁻¹
    ħ = 1.054571817e-34    # reduced Planck constant / J s
    rHH = 1.813e-10        # methyl H–H distance / m
    P2sq = 0.25            # [P₂(cos 90°)]² for an ideal methyl
    (9 / 10) * P2sq * μ0_4π^2 * γH^4 * ħ^2 / rHH^6
end

"Convert a cross-correlated relaxation rate η (s⁻¹) to S²τc in ns (eq 1)."
eta_to_s2tc_ns(η) = η / METHYL_K * 1e9

"""
    MethylCCRModel <: ParametricModel

Ratio model (eq 7) for methyl CCR buildup/decay analysis. Named to disambiguate
from the single-delay [`CCRExperiment`](@ref)/`ccr2d` analysis. Carries the fixed
prefactor `C` and the N relaxation delays `times`; `func(T, p)` evaluates eq 7 with
`p = [η, δ]`.
"""
struct MethylCCRModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    C::Float64
    times::Vector{Float64}
end

function MethylCCRModel(C, times)
    func = (T, p) -> begin
        s = sqrt(p[1]^2 + p[2]^2)
        @. C * p[1] * tanh(s * T) / (s - p[2] * tanh(s * T))
    end
    return MethylCCRModel(func, ["eta", "delta"], "Relaxation delay T / s",
                          Float64(C), collect(Float64, times))
end

"""
    methylccr2d(buildupexpt, decayexpt, T; C=3/4, skipplanes=nothing)

Start an interactive GUI for methyl ¹H–¹H cross-correlated relaxation analysis.

For each peak, the ratio of buildup to decay intensities |Iₐ/I_b| is measured across a
series of relaxation delays `T` and fitted to eq 7 with two parameters: the
cross-correlated relaxation rate `η` (s⁻¹) and a coupling term `δ` (< 0). The fitted
`η` is converted to the methyl order parameter × tumbling time `S²τc` (ns) via eq 1
(ideal methyl geometry), which is the parameter shown in the summary plot.

# Arguments
- `buildupexpt`: buildup series (Iₐ) — a single pseudo-3D path string, or a
  `Vector{String}` of per-delay 2D data directories.
- `decayexpt`: decay series (I_b), in the same form as `buildupexpt`.
- `T`: vector of relaxation delays in **seconds**, or a path string to a text file of
  delays (one per line; lines beginning with `#` are ignored). Each series must have one
  plane per delay.

# Keyword Arguments
- `C`: fixed prefactor in eq 7. `3/4` (default) for triple-quantum (TQ); `1/2` for
  double-quantum (DQ).
- `skipplanes`: optional list of delay indices (1-based, into `T`) to exclude from the
  eq 7 fit. Skipped points appear as open grey markers.

# Example
```julia
# buildup and decay each as a pseudo-3D dataset
methylccr2d("11/pdata/1", "12/pdata/1", [0.001, 0.002, 0.004, 0.006, 0.010])

# double-quantum variant
methylccr2d("11/pdata/1", "12/pdata/1", "vdlist.txt"; C=1/2)
```
"""
function methylccr2d(buildupexpt, decayexpt, T; C=3 / 4, skipplanes=nothing)
    # Parse the relaxation delays (vector, file path, or scalar)
    tau = Float64[]
    if T isa AbstractString
        append!(tau, vec(readdlm(T; comments=true)))
    elseif T isa AbstractVector
        for t in T
            if t isa AbstractString
                append!(tau, vec(readdlm(t; comments=true)))
            else
                append!(tau, t)
            end
        end
    else
        push!(tau, Float64(T))
    end
    N = length(tau)
    N > 0 || error("No relaxation delays provided")

    skip = isnothing(skipplanes) ? Int[] : collect(Int, skipplanes)
    if !isempty(skip)
        bad = filter(i -> i < 1 || i > N, skip)
        isempty(bad) ||
            error("skipplanes indices out of range (got $bad for $N delays)")
    end

    specdata = preparespecdata_methylccr(buildupexpt, decayexpt, N)
    peaks = Observable(Vector{Peak}())

    model = MethylCCRModel(C, tau)
    # x = [T; T] keeps the generic per-plane bookkeeping consistent; the CCR model
    # only ever uses its own `times`.
    expt = IntensityExperiment(specdata, peaks, model, [tau; tau],
                               ModelFitVisualisation(); skipplanes=skip)

    return gui!(expt)
end

# Load one buildup/decay series (pseudo-3D path or vector of 2D paths) into flat
# per-plane vectors, reusing the intensity-experiment loader.
function _load_ccr_series(input)
    specs = Any[]
    xs = Any[]
    ys = Any[]
    zs = Any[]
    σs = Float64[]
    files = input isa AbstractString ? [input] : collect(input)
    for f in files
        spec, x, y, z, σ = loadspecdata(f, IntensityExperiment)
        for zi in z
            push!(specs, spec)
            push!(xs, x)
            push!(ys, y)
            push!(zs, collect(zi))
            push!(σs, σ)
        end
    end
    return specs, xs, ys, zs, σs
end

# Build a 2N-plane SpecData from buildup (planes 1:N) and decay (planes N+1:2N). All
# planes are divided by a single common σ (first buildup plane) so it cancels in the
# Iₐ/I_b ratio — mirroring preparespecdata(.., CCRExperiment).
function preparespecdata_methylccr(buildupexpt, decayexpt, ntimes)
    bs = _load_ccr_series(buildupexpt)
    ds = _load_ccr_series(decayexpt)

    length(bs[4]) == ntimes ||
        error("buildup series has $(length(bs[4])) planes but $ntimes delays were given")
    length(ds[4]) == ntimes ||
        error("decay series has $(length(ds[4])) planes but $ntimes delays were given")

    specs = [bs[1]; ds[1]]
    x = [bs[2]; ds[2]]
    y = [bs[3]; ds[3]]
    z = [bs[4]; ds[4]]
    σ = [bs[5]; ds[5]]
    zlabels = [["buildup $i" for i in 1:ntimes]; ["decay $i" for i in 1:ntimes]]

    σ1 = σ[1]
    return SpecData(specs, x, y, z ./ σ1, σ ./ σ1, zlabels)
end

# Post-parameters: fitted η, δ and the derived S²τc.
function setup_post_parameters!(peak::Peak, ::MethylCCRModel)
    peak.postparameters[:eta] = Parameter("eta", 0.0)
    peak.postparameters[:delta] = Parameter("delta", 0.0)
    peak.postparameters[:S2tc] = Parameter("S2tc", 0.0)
    return
end

# Split the 2N amplitude vector into buildup (Iₐ) / decay (I_b) and return the ratio
# |Iₐ/I_b| per delay as `Measurement`s, propagating the peak-intensity uncertainties via
# Measurements.jl (the ratio error grows as the intensities become small).
function _ccr_ratio(peak, model::MethylCCRModel)
    N = length(model.times)
    amp = peak.parameters[:amp].value[]
    err = peak.parameters[:amp].uncertainty[]
    A = amp .± err
    Ia = A[1:N]
    Ib = A[(N + 1):(2N)]
    return abs.(Ia ./ Ib)
end

# Initial guess for [η, δ]. ratio ≈ C·η·T at small T, so estimate η from the
# earliest delay; δ starts at a small negative value (paper: δ < 0).
function _ccr_initial(times, ratio, C)
    pos = findall(>(0), times)
    if isempty(pos)
        η0 = 100.0
    else
        j = pos[argmin(times[pos])]
        η0 = abs(ratio[j]) / (C * times[j])
    end
    η0 = clamp(η0, 1.0, 1e5)
    return [η0, -1.0]
end

function postfit!(peak::Peak, expt::IntensityExperiment, model::MethylCCRModel)
    @debug "Post-fitting methyl CCR model"
    times = model.times
    N = length(times)
    ratio = _ccr_ratio(peak, model)

    skip = Set(expt.skipplanes)
    keep = [i for i in 1:N if i ∉ skip]
    t = times[keep]
    rval = Measurements.value.(ratio[keep])
    rerr = Measurements.uncertainty.(ratio[keep])

    # Weighted fit (wt = 1/σ²): down-weight low-intensity points whose ratio is poorly
    # determined. eps() guards against a zero uncertainty.
    wt = 1 ./ max.(rerr, eps()) .^ 2
    p0 = _ccr_initial(t, rval, model.C)
    fit = curve_fit(model.func, t, rval, wt, p0; lower=[0.0, -Inf], upper=[Inf, 0.0])
    ηm = coef(fit)[1] ± stderror(fit)[1]
    δm = coef(fit)[2] ± stderror(fit)[2]
    s2tcm = eta_to_s2tc_ns(ηm)  # eq 1, with the η uncertainty propagated through

    peak.postparameters[:eta].value[] .= Measurements.value(ηm)
    peak.postparameters[:eta].uncertainty[] .= Measurements.uncertainty(ηm)
    peak.postparameters[:delta].value[] .= Measurements.value(δm)
    peak.postparameters[:delta].uncertainty[] .= Measurements.uncertainty(δm)
    peak.postparameters[:S2tc].value[] .= Measurements.value(s2tcm)
    peak.postparameters[:S2tc].uncertainty[] .= Measurements.uncertainty(s2tcm)

    peak.postfitted[] = true
    return
end

function get_model_data(peak, expt::IntensityExperiment, model::MethylCCRModel)
    isnothing(peak) &&
        return (Point2f[], _empty_errorbars(), Point2f[], Point2f[], _empty_errorbars())

    times = model.times
    N = length(times)
    ratio = _ccr_ratio(peak, model)
    rval = Measurements.value.(ratio)
    rerr = Measurements.uncertainty.(ratio)

    skip = _skipset(expt)
    active = [i for i in 1:N if i ∉ skip]
    skipped = [i for i in 1:N if i ∈ skip]

    obs_points = Point2f.(times[active], rval[active])
    obs_errors = [(times[i], rval[i], rerr[i]) for i in active]
    skip_points = Point2f.(times[skipped], rval[skipped])
    skip_errors = [(times[i], rval[i], rerr[i]) for i in skipped]

    if peak.postfitted[]
        η = peak.postparameters[:eta].value[][1]
        δ = peak.postparameters[:delta].value[][1]
        xpred = collect(range(0.0, 1.1 * maximum(times), 100))
        ypred = model.func(xpred, [η, δ])
        fit_points = Point2f.(xpred, ypred)
    else
        fit_points = Point2f[]
    end

    return (obs_points, obs_errors, fit_points, skip_points, skip_errors)
end

function model_parameter_text(peak::Peak, ::MethylCCRModel)
    η = peak.postparameters[:eta]
    δ = peak.postparameters[:delta]
    s2tc = peak.postparameters[:S2tc]
    return ["η: $(η.value[][1] ± η.uncertainty[][1]) s⁻¹",
            "δ: $(δ.value[][1] ± δ.uncertainty[][1]) s⁻¹",
            "S²τc: $(s2tc.value[][1] ± s2tc.uncertainty[][1]) ns"]
end
