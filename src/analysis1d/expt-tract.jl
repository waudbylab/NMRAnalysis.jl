# TRACT: a TROSY / anti-TROSY pair, yielding the rotational correlation time τc.
#
#   1. entry point   2. type   3. interface   4. science   5. presentation

# ---- 1. entry point -----------------------------------------------------------

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

    B0 = 2π * acqus(trosy, :bf1) / GAMMA_H
    ωN = 2π * acqus(trosy, :bf3)
    f = tract_f(; B0)

    ds = Dataset1D(Planes(traces, vars), default_noise_center(trosy))
    expt = isnothing(regions) ? TractExperiment(ds; ωN, f) :
           TractExperiment(ds; ωN, f, regions)
    return run1d(expt; integration)
end

# ---- 2. type ------------------------------------------------------------------

"""
    TractExperiment(dataset; ωN, f, regions=…)

Fit TROSY and anti-TROSY decays (grouped by `vars.which ∈ {:trosy, :anti}`) and derive
the rotational correlation time τc from the cross-correlated relaxation-rate difference
ΔR = R(anti) − R(trosy). `ωN` is the ¹⁵N Larmor frequency (rad s⁻¹) and `f` the
dipole/CSA cross-correlation prefactor (see `tract_f`).
"""
struct TractExperiment <: Experiment1D
    dataset::Dataset1D
    regions::Vector{Region}
    ωN::Float64
    f::Float64
end

function TractExperiment(dataset::Dataset1D; ωN, f, regions=[defaultregion(dataset)])
    return TractExperiment(dataset, collect(Region, regions), Float64(ωN), Float64(f))
end

# ---- 3. interface -------------------------------------------------------------

seriesmodel(::TractExperiment) = ExponentialModel()
fitaxis(::TractExperiment) = :time
groupcols(::TractExperiment) = (:which,)
primaryparam(::TractExperiment) = :τc

# ---- 4. science ---------------------------------------------------------------

"¹H gyromagnetic ratio / rad s⁻¹ T⁻¹."
const GAMMA_H = 2.6752218744e8

function postfitglobal!(results::AbstractVector{RegionResult}, e::TractExperiment)
    # Iterate the region labels actually present in `results`, not `regions(e)` (the
    # experiment's own, fixed-at-construction region list) - the GUI's live `regs`
    # argument to `analyse` can include regions added interactively after construction,
    # and those still need a τc summary.
    for label in unique(r.region for r in results)
        rs = filter(r -> r.region == label, results)
        trosy = findfirst(r -> r.group.which == :trosy, rs)
        anti = findfirst(r -> r.group.which == :anti, rs)
        (isnothing(trosy) || isnothing(anti)) && continue
        ηxy = (param(rs[anti], :R) - param(rs[trosy], :R)) / 2
        τc = tract_tauc(e.f, e.ωN, ηxy)
        # Recorded on both series of the pair: either is a complete answer for the
        # region, and the GUI shows whichever the user has selected.
        for r in (rs[trosy], rs[anti])
            setpost!(r, :ηxy, ηxy)
            setpost!(r, :τc, τc)
        end
    end
    return nothing
end

"""
    tract_f(; B0, θ=17π/180) -> Float64

Dipole/CSA cross-correlation prefactor used in the TRACT τc relation, given the static
field `B0` (T). Constants follow the standard ¹⁵N–¹H amide treatment.
"""
function tract_f(; B0, θ=17 * π / 180)
    μ0 = 4π * 1e-7
    γH = GAMMA_H
    γN = -2.7126180e7
    ħ = 6.62607015e-34 / 2π
    rNH = 1.02e-10
    ΔδN = 160e-6
    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = B0 * γN * ΔδN / (3 * sqrt(2))
    return p * c * (3cos(θ)^2 - 1)
end

"""
    tract_tauc(f, ωN, ηxy) -> Float64

Rotational correlation time τc (ns) from the cross-correlated cross-relaxation rate
`ηxy`, by the analytic inversion of `ηxy = f·(4/5·τc + 3/5·τc/(1+(ωN·τc)²))` used in the
existing `tract` routine.
"""
function tract_tauc(f, ωN, ηxy)
    x = sqrt(21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy^2 * ωN^8 + 625 * f^2 * ηxy^4 * ωN^10)
    y = cbrt(1800 * f^2 * ηxy * ωN^4 + 125 * ηxy^3 * ωN^6 + 24 * sqrt(3) * x)
    τc = (5 * ηxy) / (24 * f) -
         (336 * f^2 * ωN^2 - 25 * ηxy^2 * ωN^4) / (24 * f * ωN^2 * y) + y / (24 * f * ωN^2)
    return 1e9 * τc
end

# ---- 5. presentation ----------------------------------------------------------

windowtitle(::TractExperiment) = "TRACT"

resultxfactor(e::TractExperiment) = timescale(column(dataset(e).planes, :time))[1]

function result_labels(e::TractExperiment)
    _, unit = timescale(column(dataset(e).planes, :time))
    return ("Relaxation delay / $unit", "Integrated intensity (a.u.)")
end

seriesnames(::TractExperiment) = ["TROSY", "anti-TROSY"]

function spectruminfo(::TractExperiment, vars::NamedTuple)
    which = vars.which == :trosy ? "TROSY" : "anti-TROSY"
    return "$(round(vars.time; digits=3)) s delay ($which)"
end
