# B₁ fields: how a power level maps to a field strength, and how that field is distributed
# across the sample. Everything about B₁ lives here, so the whole story can be read in one
# sitting.
#
# This file is deliberately free of any dependency beyond NMRTools and Measurements, so it
# can be promoted to NMRTools (where `Power`, `db`, `hz` and `referencepulse` already live)
# once the design has settled. Fitting nutation *data* is not here: that belongs to the
# analysis, in `analysis1d/expt-nutation.jl`, which builds a `B1Calibration` from its fits.

"""
    DEFAULTINHOMOGENEITY

Fractional B₁ inhomogeneity assumed when nothing has been measured: 5%, typical of a
standard probe. Only a nutation calibration can do better than this guess.
"""
const DEFAULTINHOMOGENEITY = 0.05

"""
    B1Calibration(power, ν1; nuc=nothing, inhomogeneity=DEFAULTINHOMOGENEITY, source="")
    B1Calibration(spec, nuc; inhomogeneity=DEFAULTINHOMOGENEITY)

The B₁ field a power level produces, and how widely that field is distributed across the
sample.

The field follows the ideal square-root power law,

    ν₁(p) = ν₁ref · 10^(−L·(dB(p) − dB(pref))/20),

in which the `linearity` L is 1 for a perfectly linear amplifier. One measurement fixes
`ν₁ref` and assumes L = 1; two or more fit both, so L is a direct diagnostic of amplifier
compression, and departures from 1 are what a multi-power calibration is for.

The second form is the *nominal* calibration read from a spectrum's own reference pulse
(`p1`/`pl1`, via `referencepulse`): one point, L = 1, and an assumed inhomogeneity. It
reproduces `hz(p, refpower, refpulse, 90)` exactly, so an analysis given no calibration
behaves as it always did.

Use `hz(power, cal)` to get a field from a power and `Power(ν1, cal)` to invert it.

# Arguments
- `power`: the power levels calibrated at, as NMRTools `Power`s.
- `ν1`: the field measured at each, in Hz, with uncertainties.
- `nuc`: the nucleus (channel) this calibration applies to, for checking it is applied to
  the right one.
- `inhomogeneity`: fractional width of the B₁ distribution (0.05 for 5%).
- `source`: where the calibration came from, for provenance in saved results.

# Example
```julia
cal = B1Calibration([Power(18.0, :dB), Power(12.0, :dB)], [1250.0 ± 5.0, 2480.0 ± 9.0])
hz(Power(15.0, :dB), cal)    # field at an intermediate power, Hz
Power(2000.0, cal)           # power needed for a 2 kHz field
```
"""
struct B1Calibration
    nuc::Any
    power::Vector{Power{Float64}}
    ν1::Vector{Measurement{Float64}}
    powerref::Power{Float64}
    ν1ref::Measurement{Float64}
    linearity::Measurement{Float64}
    inhomogeneity::Float64
    source::String
end

function B1Calibration(power::AbstractVector, ν1::AbstractVector;
                       nuc=nothing, inhomogeneity::Real=DEFAULTINHOMOGENEITY,
                       source::AbstractString="")
    length(power) == length(ν1) ||
        throw(ArgumentError("got $(length(power)) powers for $(length(ν1)) field strengths"))
    isempty(power) && throw(ArgumentError("a calibration needs at least one measurement"))
    0 ≤ inhomogeneity < 1 ||
        throw(ArgumentError("inhomogeneity is a fraction, not a percentage " *
                            "(got $inhomogeneity)"))
    pow = collect(Power{Float64}, power)
    ν = [x isa Measurement ? x : x ± 0.0 for x in ν1]
    all(>(0), Measurements.value.(ν)) ||
        throw(ArgumentError("field strengths must be positive: got " *
                            "$(Measurements.value.(ν)) Hz"))
    # Anchored on the first measurement, so a single-point calibration carries the measured
    # value itself rather than one round-tripped through a logarithm.
    νref, lin = length(pow) == 1 ? (first(ν), 1.0 ± 0.0) :
                calibrationfit(db.(pow), ν, db(first(pow)))
    return B1Calibration(nuc, pow, ν, first(pow), νref, lin,
                         Float64(inhomogeneity), String(source))
end

function B1Calibration(spec, nuc; inhomogeneity::Real=DEFAULTINHOMOGENEITY)
    pulse, pow = referencepulse(spec, nuc)
    # 90 / (360 * pulse), written as NMRTools' own `hz` writes it, so that the nominal
    # calibration and the four-argument `hz` agree to the last bit.
    return B1Calibration([pow], [90 / (360 * pulse)]; nuc, inhomogeneity,
                         source="reference pulse")
end

"""
    calibrationfit(dB, ν1, dBref) -> (ν1ref, linearity)

Weighted straight-line fit of log₁₀ν₁ against power, returning the field at `dBref` and the
linearity (1 for an ideal amplifier).

The fit is linear because the ideal power law is a straight line in these coordinates, so
there is nothing to iterate and no initial guess to get wrong. Written as the five weighted
sums of the standard treatment (Numerical Recipes §15.2), with the abscissa centred on
`dBref` so that the intercept *is* log₁₀ν₁ref and its variance comes straight from the
covariance matrix. Uncertainties on ν₁ transform as σ(log₁₀ν) = σ(ν)/(ν ln10).
"""
function calibrationfit(dB::AbstractVector, ν1::AbstractVector, dBref::Real)
    allequal(dB) &&
        throw(ArgumentError("cannot fit a calibration curve: every measurement is at the " *
                            "same power ($(first(dB)) dB)"))
    ν = Measurements.value.(ν1)
    all(>(0), ν) || throw(ArgumentError("field strengths must be positive: got $ν Hz"))
    y = log10.(ν)
    σ = Measurements.uncertainty.(ν1) ./ (ν .* log(10))
    w = all(>(0), σ) ? 1 ./ σ .^ 2 : ones(length(y))
    x = dB .- dBref

    S = sum(w)
    Sx = sum(w .* x)
    Sxx = sum(w .* x .^ 2)
    Sy = sum(w .* y)
    Sxy = sum(w .* x .* y)
    Δ = S * Sxx - Sx^2
    c0 = (Sxx * Sy - Sx * Sxy) / Δ
    c1 = (S * Sxy - Sx * Sy) / Δ

    # Standard errors from the covariance matrix; scaled by the reduced chi-squared where
    # the points are unweighted, there being no other estimate of their scatter then.
    scale = all(>(0), σ) ? 1.0 : sum(w .* (y .- c0 .- c1 .* x) .^ 2) / max(length(y) - 2, 1)
    ν0 = 10^c0
    return (ν0 ± ν0 * log(10) * sqrt(scale * Sxx / Δ), -20c1 ± 20sqrt(scale * S / Δ))
end

"""
    hz(power, calibration) -> Float64

The B₁ field strength (Hz) a `power` level produces, per `calibration`.

The calibration's own uncertainty is not propagated: it is far smaller than the B₁
inhomogeneity that [`B1Distribution`](@ref) accounts for, and the analyses downstream treat
the nominal field as exact.
"""
function NMRTools.hz(p::Power, cal::B1Calibration)
    ΔdB = Measurements.value(linearity(cal)) * (db(p) - db(refpower(cal)))
    return Measurements.value(ν1ref(cal)) * 10^(-ΔdB / 20)
end

"""
    Power(ν1, calibration) -> Power

The power level giving a field strength of `ν1` Hz, inverting `hz(power, cal)`. This is
what turns a list of wanted spin-lock strengths into a list of powers to set.
"""
function NMRTools.Power(ν1::Real, cal::B1Calibration)
    ν1 > 0 || throw(ArgumentError("field strength must be positive (got $ν1 Hz)"))
    ratio = ν1 / Measurements.value(ν1ref(cal))
    ΔdB = -20 * log10(ratio) / Measurements.value(linearity(cal))
    return Power(db(refpower(cal)) + ΔdB, :dB)
end

# A calibration is one object, not a container of its measurements, so `hz.(powers, cal)`
# broadcasts over the powers alone.
Base.broadcastable(cal::B1Calibration) = Ref(cal)

"The power levels a calibration was measured at."
powers(cal::B1Calibration) = cal.power

"The field strength measured at each power level, with its uncertainty."
fields(cal::B1Calibration) = cal.ν1

"Field strength (Hz) at the calibration's reference power, with its uncertainty."
ν1ref(cal::B1Calibration) = cal.ν1ref

"The calibration's reference power."
refpower(cal::B1Calibration) = cal.powerref

"""
    linearity(calibration)

Fitted power-law exponent relative to the ideal, so 1 means the amplifier follows
ν₁ ∝ √W. Below 1 the field falls short of the ideal law as the power is raised, the usual
signature of amplifier compression. Exactly 1, without uncertainty, where a single
measurement left nothing to fit.
"""
linearity(cal::B1Calibration) = cal.linearity

"""
    inhomogeneity(calibration) -> Float64

Fractional width of the B₁ distribution (0.05 for 5%).
"""
inhomogeneity(cal::B1Calibration) = cal.inhomogeneity

"The nucleus (channel) a calibration was measured on, or `nothing` if not recorded."
NMRTools.nucleus(cal::B1Calibration) = cal.nuc

"Where the calibration came from, for provenance in saved results."
source(cal::B1Calibration) = cal.source

npoints(cal::B1Calibration) = length(cal.power)

function Base.show(io::IO, cal::B1Calibration)
    print(io, "B1Calibration(")
    isnothing(nucleus(cal)) || print(io, nucleus(cal), ", ")
    print(io, npoints(cal), npoints(cal) == 1 ? " point" : " points", ", ")
    print(io, round(Measurements.value(ν1ref(cal)); sigdigits=5), " Hz at ",
          round(db(refpower(cal)); digits=2), " dB")
    npoints(cal) == 1 ||
        print(io, ", linearity ", round(Measurements.value(linearity(cal)); digits=4))
    print(io, ", B₁ inhomogeneity ", round(100 * inhomogeneity(cal); digits=2), "%")
    isempty(source(cal)) || print(io, ", from ", source(cal))
    return print(io, ")")
end

# ---- the distribution of B₁ across the sample ---------------------------------

# Three-point Gauss-Hermite quadrature, in units of the distribution's standard deviation.
# These nodes and weights integrate a Gaussian exactly up to its fifth moment, so the
# sampled set has the right mean and the right width for any σ. Node *positions* have to
# scale with σ for that to hold: fixing them at, say, 0.97/1.00/1.05 of nominal caps the
# representable width at 3.9% however large σ is said to be.
const B1NODES = (-√3.0, 0.0, √3.0)
const B1WEIGHTS = (1 / 6, 2 / 3, 1 / 6)

"""
    B1Distribution(σ)
    B1Distribution(calibration)
    B1Distribution(scaling, weight)

The distribution of B₁ across the sample, as the quadrature a simulation actually needs: a
set of field `scaling`s relative to nominal and the `weight` each carries.

`σ` is the fractional width of a Gaussian distribution, sampled at the three Gauss-Hermite
nodes `1 ± √3σ` with weights 1/6, 2/3, 1/6. `σ = 0` gives the single nominal field, so a
simulation with no inhomogeneity costs exactly what it did before.

Holding the quadrature rather than the distribution is what keeps this general: another
family, or a B₁ profile measured on the probe, is a different set of scalings and weights
and nothing downstream changes. The third form takes them directly.

# Example
```julia
d = B1Distribution(0.05)          # 91.3%, 100% and 108.7% of nominal
b1average(ν -> ν^2, d, 1000.0)    # ⟨ν²⟩ over the distribution
```
"""
struct B1Distribution
    scaling::Vector{Float64}
    weight::Vector{Float64}
    function B1Distribution(scaling, weight)
        length(scaling) == length(weight) ||
            throw(ArgumentError("got $(length(scaling)) scalings for " *
                                "$(length(weight)) weights"))
        all(>(0), scaling) ||
            throw(ArgumentError("field scalings must be positive: got $scaling"))
        total = sum(weight)
        total > 0 || throw(ArgumentError("weights must sum to a positive number"))
        return new(collect(float.(scaling)), collect(float.(weight)) ./ total)
    end
end

function B1Distribution(σ::Real)
    σ ≥ 0 || throw(ArgumentError("B₁ inhomogeneity cannot be negative (got $σ)"))
    σ < 1 / √3 ||
        throw(ArgumentError("a B₁ inhomogeneity of $σ is too wide for a three-point " *
                            "quadrature: the lowest node, 1 - √3σ, is at or below zero"))
    iszero(σ) && return B1Distribution([1.0], [1.0])
    return B1Distribution(1 .+ σ .* B1NODES, B1WEIGHTS)
end

B1Distribution(cal::B1Calibration) = B1Distribution(inhomogeneity(cal))

"""
    b1average(f, distribution, ν1)

Average `f(ν)` over the B₁ `distribution` of a nominal field `ν1`.

This is the right average for anything linear in the observable, which is to say for a
directly measured intensity: a CEST profile is the weighted sum of the profiles the
separate parts of the sample give. For a quantity that reaches the data through a fit, see
[`b1rate`](@ref).
"""
function b1average(f, d::B1Distribution, ν1)
    return sum(w * f(s * ν1) for (s, w) in zip(d.scaling, d.weight))
end

"""
    b1rate(f, distribution, ν1, T)

Apparent relaxation rate over a B₁ `distribution`, where `f(ν)` is the rate at field `ν`
and `T` the evolution time the rate was measured over.

A spread of B₁ gives a spread of rates, so the decay is a sum of exponentials, and what a
monoexponential fit to the data reports is

    R = −ln(Σᵢ wᵢ exp(−Rᵢ T)) / T

rather than the weighted mean of the rates ⟨R⟩. The two agree as T → 0 and diverge by
about ½T·var(R) beyond that, which at the low spin-lock strengths of an R₁ρ dispersion is
comparable with the uncertainty on the measurement. `T` is the time the experiment actually
sampled, so the simulated quantity is the one the data reduction produced.
"""
function b1rate(f, d::B1Distribution, ν1, T)
    T > 0 || return b1average(f, d, ν1)
    return -log(b1average(ν -> exp(-f(ν) * T), d, ν1)) / T
end
