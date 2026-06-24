using NMRAnalysis
using Test
using Statistics
import Measurements

const A1 = NMRAnalysis.Analysis1D

# ---- helpers ------------------------------------------------------------------

"Lorentzian peak shape centred at δ0."
lorentz(δ, δ0, w) = @. 1 / (1 + ((δ - δ0) / w)^2)

"A trace built from a list of (δ0, amplitude) peaks on a shared grid."
function make_trace(δ, peaks; w=0.05)
    y = zeros(length(δ))
    for (δ0, amp) in peaks
        y .+= amp .* lorentz(δ, δ0, w)
    end
    return A1.Trace(δ, y)
end

const δgrid = collect(range(0.0, 10.0; length=2000))
const noiseR = A1.Region("noise", 9.0, 9.5)   # empty region ⇒ σ = 0 ⇒ unweighted fits

# ===============================================================================
@testset "Analysis1D: core reductions" begin
    # boxcar: exactly 2.0 in [4.5,5.5], 0 elsewhere — gives a flat (zero-σ) noise window
    t = A1.Trace(δgrid, [4.5 ≤ δ ≤ 5.5 ? 2.0 : 0.0 for δ in δgrid])
    roi = A1.Region("p", 4.5, 5.5)

    @test A1.integrate(t, roi) > 0
    # a zero-width region falls back to the nearest point (a height)
    height = A1.Region("h", 5.0)
    @test A1.integrate(t, height) ≈ 2.0 atol = 1e-6

    # noise propagation: σ = 0 over the flat (all-zero) noise region
    ds = A1.Dataset1D(A1.Planes([t], [(; time=0.0)]), A1.Region("noise", 9.0, 9.5))
    I = A1.integrals(roi, ds)
    @test length(I) == 1
    @test Measurements.uncertainty(I[1]) == 0.0
end

@testset "Analysis1D: grouping" begin
    vars = [(; time=0.0, which=:trosy), (; time=1.0, which=:trosy),
            (; time=0.0, which=:anti)]
    planes = A1.Planes([make_trace(δgrid, [(5.0, 1.0)]) for _ in 1:3], vars)
    groups = A1.groupseries(planes, (:which,))
    @test length(groups) == 2
    @test A1.column(planes, :which) == [:trosy, :trosy, :anti]
end

# ===============================================================================
@testset "Analysis1D: relaxation recovers R" begin
    Rtrue = 12.0
    times = collect(range(0.0, 0.3; length=12))
    traces = [make_trace(δgrid, [(5.0, exp(-Rtrue * t))]) for t in times]
    ds = A1.Dataset1D(A1.Planes(traces, [(; time=t) for t in times]), noiseR)

    expt = A1.RelaxationExperiment(ds; regions=[A1.Region("p", 4.5, 5.5)])
    res = A1.analyse(expt)
    @test length(res.series) == 1
    R = A1.param(res.series[1], "R")
    @test Measurements.value(R) ≈ Rtrue rtol = 1e-3
end

@testset "Analysis1D: recovery (IR) recovers R" begin
    Rtrue = 1.5
    times = collect(range(0.0, 3.0; length=15))
    amp(t) = 1 - 2 * exp(-Rtrue * t)              # full inversion, A=1, C=2
    # shift baseline positive so the Lorentzian amplitude stays sensible under integration
    traces = [make_trace(δgrid, [(5.0, amp(t) + 2.0)]) for t in times]
    ds = A1.Dataset1D(A1.Planes(traces, [(; time=t) for t in times]), noiseR)

    expt = A1.RelaxationExperiment(ds; ir=true, regions=[A1.Region("p", 4.5, 5.5)])
    res = A1.analyse(expt)
    R = A1.param(res.series[1], "R")
    @test Measurements.value(R) ≈ Rtrue rtol = 1e-2
end

# ===============================================================================
@testset "Analysis1D: TRACT τc" begin
    # tract_tauc monotonic in ηxy, and positive
    f = A1.tract_f(; B0=14.1)
    ωN = 2π * 60.8e6
    @test A1.tract_tauc(f, ωN, 5.0) > 0
    @test A1.tract_tauc(f, ωN, 10.0) > A1.tract_tauc(f, ωN, 5.0)

    # full pipeline: anti-TROSY relaxes faster than TROSY
    times = collect(range(0.0, 0.05; length=10))
    Rt, Ra = 20.0, 45.0
    tr = [make_trace(δgrid, [(8.0, exp(-Rt * t))]) for t in times]
    at = [make_trace(δgrid, [(8.0, exp(-Ra * t))]) for t in times]
    vars = vcat([(; time=t, which=:trosy) for t in times],
                [(; time=t, which=:anti) for t in times])
    ds = A1.Dataset1D(A1.Planes(vcat(tr, at), vars), noiseR)

    expt = A1.TractExperiment(ds; ωN=ωN, f=f, regions=[A1.Region("p", 7.5, 8.5)])
    res = A1.analyse(expt)
    s = res.summary[1]
    @test Measurements.value(s.Ranti) > Measurements.value(s.Rtrosy)
    @test Measurements.value(s.τc) > 0
end

# ===============================================================================
@testset "Analysis1D: nutation 90° pulse" begin
    νtrue = 250.0                                  # Hz
    durs = collect(range(0.0, 0.002; length=25))   # estimate 0.5/range ≈ νtrue
    traces = [make_trace(δgrid, [(5.0, sin(2π * νtrue * t) * exp(-30 * t) + 1e-9)])
              for t in durs]
    ds = A1.Dataset1D(A1.Planes(traces, [(; duration=d) for d in durs]), noiseR)

    expt = A1.NutationExperiment(ds; regions=[A1.Region("p", 4.8, 5.2)])
    res = A1.analyse(expt)
    ν = res.summary[1].ν
    @test Measurements.value(ν) ≈ νtrue rtol = 1e-2
    @test Measurements.value(res.summary[1].pulse90) ≈ 1 / (4νtrue) rtol = 1e-2
end

# ===============================================================================
@testset "Analysis1D: STD contrast, buildup, epitope" begin
    # --- single saturation time: STD% = (Iref - Ion)/Iref ---
    L = A1.Region("L1", 2.8, 3.2)
    ref = make_trace(δgrid, [(3.0, 1.0)])
    on = make_trace(δgrid, [(3.0, 0.8)])           # 20% reduction
    ds = A1.Dataset1D(A1.Planes([ref, on],
                                [(; sat=:reference, tsat=0.5), (; sat=:methyl, tsat=0.5)]),
                      noiseR)
    expt = A1.STDExperiment(ds; regions=[L])
    out = A1.analyse(expt)
    @test length(out.points) == 1
    @test Measurements.value(out.points[1].std) ≈ 0.2 rtol = 1e-3

    # --- buildup: STD(t) = STDmax (1 - exp(-k t)) ---
    STDmax, k = 0.3, 5.0
    tsats = [0.1, 0.3, 0.6, 1.0]
    traces = A1.Trace[]
    vars = NamedTuple[]
    for τ in tsats
        push!(traces, make_trace(δgrid, [(3.0, 1.0)]));              push!(vars, (; sat=:reference, tsat=τ))
        push!(traces, make_trace(δgrid, [(3.0, 1 - STDmax * (1 - exp(-k * τ)))])); push!(vars, (; sat=:methyl, tsat=τ))
    end
    ds2 = A1.Dataset1D(A1.Planes(traces, vars), noiseR)
    out2 = A1.analyse(A1.STDExperiment(ds2; regions=[L]))
    @test length(out2.buildups) == 1
    b = out2.buildups[1]
    @test Measurements.value(b.std_af_max) ≈ STDmax rtol = 1e-2
    @test Measurements.value(b.k) ≈ k rtol = 5e-2
    @test Measurements.value(b.std_af0) ≈ STDmax * k rtol = 5e-2

    # --- epitope normalisation across two regions ---
    L1 = A1.Region("L1", 2.8, 3.2)
    L2 = A1.Region("L2", 5.8, 6.2)
    reft = make_trace(δgrid, [(3.0, 1.0), (6.0, 1.0)])
    ont = make_trace(δgrid, [(3.0, 0.7), (6.0, 0.85)])    # STD 0.30 and 0.15
    ds3 = A1.Dataset1D(A1.Planes([reft, ont],
                                 [(; sat=:reference, tsat=0.5), (; sat=:methyl, tsat=0.5)]),
                       noiseR)
    out3 = A1.analyse(A1.STDExperiment(ds3; regions=[L1, L2]))
    rel = Dict(e.region => e.relative for e in out3.epitope)
    @test rel["L1"] ≈ 1.0 atol = 1e-6
    @test rel["L2"] ≈ 0.5 rtol = 5e-2
end

# ===============================================================================
@testset "Analysis1D: kinetics intensity vs time" begin
    times = [0.0, 1.0, 2.0, 3.0]
    traces = [make_trace(δgrid, [(4.0, exp(-0.5t))]) for t in times]

    # single run
    ds = A1.Dataset1D(A1.Planes(traces, [(; time=t) for t in times]), noiseR)
    expt = A1.KineticsExperiment(ds; regions=[A1.Region("p", 3.5, 4.5)])
    res = A1.analyse(expt)
    @test length(res.series) == 1
    @test length(res.series[1].y) == 4

    # two runs ⇒ grouped into two series
    traces2 = vcat(traces, traces)
    vars2 = vcat([(; time=t, run=1) for t in times], [(; time=t, run=2) for t in times])
    ds2 = A1.Dataset1D(A1.Planes(traces2, vars2), noiseR)
    res2 = A1.analyse(A1.KineticsExperiment(ds2; regions=[A1.Region("p", 3.5, 4.5)]))
    @test length(res2.series) == 2
end
