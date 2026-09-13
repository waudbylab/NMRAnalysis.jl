# Analysis1D: the science, with no GUI and no data files.
#
# Every test here builds its spectra in memory as `Trace`s, which is the seam the module's
# "keep the science pure" split exists to provide: `Dataset1D` → `analyse` needs no
# NMRData, no Bruker directory and no region picking. Tests that need real acquisition
# parameters (annotations, vdlists, p30/d20/gpnam, plane flattening) belong in a separate
# file alongside example data, and go through the entry points with `prompt=false` and an
# `integration` triple.

using NMRAnalysis
using NMRAnalysis.Analysis1D: Trace, Planes, groupseries, hasvar, column, nplanes
using NMRAnalysis.Analysis1D: integrate
using NMRAnalysis.Analysis1D: width, centre, recentre, setwidth, defaultregionwidth,
                              defaultamideregion
using NMRAnalysis.Analysis1D: fitaxis, groupcols, primaryparam, isconverged
using NMRAnalysis.Analysis1D: relaxationmodel, relaxationseriesmodel, nutationphase,
                              shapefactor, solventname, tractf, tracttauc
using NMRAnalysis.Analysis1D: resultstable, seriestable, experimentinfo,
                              writeresults!, csvcolumn, safename, sources, seriesname,
                              baseparam, regionlisttable, readregions!, NOISE_LABEL
using NMRAnalysis.Analysis1D: analysiscall, callstring, callvalue, writesummary
using NMRAnalysis.Analysis1D: ask, askvector, askchoice, parsevector, acqusvalue
using Measurements
using Random
using Test

# ---- synthetic spectra --------------------------------------------------------

const δ = collect(range(-2.0, 12.0; length=1401))   # 0.01 ppm/point

"""Intensities of a Lorentzian peak of amplitude `a` at `δ0`, over the shared axis."""
lorentzian(a, δ0; hw=0.03) = @. a * hw^2 / (hw^2 + (δ - δ0)^2)

"""A single-peak spectrum, with deterministic noise when an `rng` is given."""
function peaktrace(a, δ0; hw=0.03, noise=0.0, rng=nothing)
    y = lorentzian(a, δ0; hw)
    isnothing(rng) || (y .+= noise .* randn(rng, length(δ)))
    return Trace(δ, y)
end

"""A dataset of one peak whose amplitude follows `amps`, tagged by `var => vals`."""
function peakdataset(amps, var::Symbol, vals; δ0=8.0, noise=0.0, seed=42, noisecentre=0.0)
    rng = MersenneTwister(seed)
    traces = [peaktrace(a, δ0; noise, rng) for a in amps]
    vars = [NamedTuple{(var,)}((Float64(v),)) for v in vals]
    return Dataset1D(Planes(traces, vars), noisecentre, "synthetic")
end

signalregion(δ0=8.0) = [Region("signal", δ0 - 0.3, δ0 + 0.3)]

@testset "Analysis1D" begin
    @testset "Region" begin
        r = Region("signal", 9.0, 7.0)
        @test r.lo == 7.0 && r.hi == 9.0        # bounds stored sorted
        @test width(r) == 2.0
        @test centre(r) == 8.0
        @test width(recentre(r, 3.0)) == 2.0
        @test centre(recentre(r, 3.0)) == 3.0
        @test width(setwidth(r, 0.5)) == 0.5
        @test centre(setwidth(r, 0.5)) == 8.0
        @test width(Region("height", 4.0)) == 0.0
        @test defaultregionwidth(δ) ≈ 0.02 * 14.0
        @test defaultamideregion().lo == 7.5
    end

    @testset "Planes and grouping" begin
        traces = [peaktrace(1.0, 8.0) for _ in 1:4]
        vars = [(; time=0.1, which=:trosy), (; time=0.2, which=:trosy),
                (; time=0.1, which=:anti), (; time=0.2, which=:anti)]
        p = Planes(traces, vars)
        @test nplanes(p) == 4
        @test column(p, :time) == [0.1, 0.2, 0.1, 0.2]
        @test hasvar(p, :which)
        @test !hasvar(p, :gradient)

        @test length(groupseries(p, ())) == 1            # ungrouped: one series
        @test last(only(groupseries(p, ()))) == 1:4

        groups = groupseries(p, (:which,))
        @test length(groups) == 2
        @test first(groups[1]) == (; which=:trosy)
        @test last(groups[1]) == [1, 2]
        @test last(groups[2]) == [3, 4]

        @test_throws ArgumentError Planes(traces, vars[1:3])
        @test_throws ArgumentError Planes(traces, [(; a=1), (; b=2), (; a=3), (; a=4)])
        @test_throws ArgumentError Trace([1.0, 2.0], [1.0])
    end

    @testset "Integration and noise" begin
        ds = peakdataset([100.0, 50.0], :time, [0.0, 1.0]; noise=1.0, noisecentre=0.0)
        I = integrate(only(signalregion()), ds)
        @test length(I) == 2
        @test Measurements.value(I[1]) > Measurements.value(I[2]) > 0
        @test Measurements.value(I[1]) ≈ 2 * Measurements.value(I[2]) rtol = 0.05
        @test Measurements.uncertainty(I[1]) > 0
        # the same σ applies to every plane: it is one measurement of the noise region
        @test Measurements.uncertainty(I[1]) == Measurements.uncertainty(I[2])

        # a zero-width region is a height: the single nearest point
        t = peaktrace(100.0, 8.0)
        @test integrate(Region("height", 8.0), t) ≈ 100.0 rtol = 1e-3
        @test integrate(Region("wide", 7.7, 8.3), t) > integrate(Region("height", 8.0), t)
    end

    @testset "Relaxation" begin
        times = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2]
        R = 2.5
        ds = peakdataset([100 * exp(-R * t) for t in times], :time, times; noise=0.5)
        expt = RelaxationExperiment(ds; regions=signalregion())
        @test fitaxis(expt) == :time
        @test primaryparam(expt) == :R

        res = analyse1d(expt)
        @test length(res) == 1
        @test res[1].region == "signal"
        @test isconverged(res[1])
        @test Measurements.value(param(res[1], :R)) ≈ R rtol = 0.05
        @test Measurements.uncertainty(param(res[1], :R)) > 0
        @test !res[1].postfitted                 # relaxation derives nothing

        # inversion recovery
        C = 2.0
        amps = [100 * (1 - C * exp(-R * t)) for t in times]
        dsr = peakdataset(amps, :time, times; noise=0.5)
        resr = analyse1d(RelaxationExperiment(dsr; model=:recovery, regions=signalregion()))
        @test Measurements.value(param(resr[1], :R)) ≈ R rtol = 0.1
        @test Measurements.value(param(resr[1], :C)) ≈ C rtol = 0.1

        @test relaxationmodel(nothing) === nothing
        @test relaxationmodel("inversion_recovery") === :recovery
        @test relaxationmodel("exponential_decay") === :exponential
        @test relaxationseriesmodel(:exponential) isa CurveFitModel
        @test relaxationseriesmodel(NoFitting()) isa NoFitting
        @test_throws ArgumentError relaxationseriesmodel(:nonsense)
    end

    @testset "TRACT" begin
        # the τc inversion, checked against the forward relation it inverts
        f = tractf(; B0=18.8)
        ωN = 2π * 81.08e6
        # ηxy = f·[4J(0) + 3J(ωN)] with J(ω) = (2/5)τc/(1+ω²τc²). Getting the (2/5)
        # wrong is a silent factor of two in every reported τc, so it is pinned here.
        for τc in (2e-9, 5e-9, 15e-9, 30e-9, 50e-9)
            ηxy = f * (8 / 5 * τc + 6 / 5 * τc / (1 + (ωN * τc)^2))
            @test tracttauc(f, ωN, ηxy) ≈ 1e9 * τc rtol = 1e-6
        end

        times = [0.0, 0.005, 0.01, 0.02, 0.04]
        Rt, Ra = 16.0, 54.0
        traces = vcat([peaktrace(100 * exp(-Rt * t), 8.0) for t in times],
                      [peaktrace(100 * exp(-Ra * t), 8.0) for t in times])
        vars = vcat([(; time=t, which=:trosy) for t in times],
                    [(; time=t, which=:anti) for t in times])
        ds = Dataset1D(Planes(traces, vars), 0.0, "synthetic")
        expt = TractExperiment(ds; ωN, f, regions=signalregion())
        @test groupcols(expt) == (:which,)

        res = analyse1d(expt)
        # one region, two series, and one set of parameters naming them apart
        @test length(res) == 1
        r = only(res)
        @test length(r.series) == 2
        @test Set(s.group.which for s in r.series) == Set([:trosy, :anti])
        @test Measurements.value(param(r, :R_trosy)) ≈ Rt rtol = 0.02
        @test Measurements.value(param(r, :R_anti)) ≈ Ra rtol = 0.02

        # η and τc are computed from both rates and sit beside them on the same region
        @test r.postfitted
        @test Measurements.value(param(r, :etaxy)) ≈ (Ra - Rt) / 2 rtol = 0.05
        @test Measurements.value(param(r, :tauc)) > 0           # ns
        @test baseparam(:R_trosy) == :R
        @test seriesname("R", (; which=:anti)) == :R_anti
        @test seriesname("R", NamedTuple()) == :R
    end

    @testset "Nutation calibration" begin
        # the initial guess assumes half a period across the sampled range, so sample one
        ν, Rdecay = 500.0, 500.0
        durations = collect(range(0.0, 1.0e-3; length=21))
        amps = [100 * sin(2π * ν * t) * exp(-Rdecay * t) for t in durations]
        ds = peakdataset(amps, :duration, durations; noise=0.2)
        expt = NutationExperiment(ds; regions=signalregion())
        @test primaryparam(expt) == :pulse90

        res = analyse1d(expt)
        @test Measurements.value(param(res[1], :nu)) ≈ ν rtol = 0.05
        # stored in the units they are quoted in: µs and %, not seconds and a fraction
        @test Measurements.value(param(res[1], :pulse90)) ≈ 1e6 / (4ν) rtol = 0.05
        @test Measurements.value(param(res[1], :inhomogeneity)) ≈
              100 * Rdecay / (2π * ν) rtol = 0.1

        @test nutationphase(nothing) === nothing
        @test nutationphase("cosine_modulated") === :cosine
        @test nutationphase("sine_modulated") === :sine
    end

    @testset "Diffusion" begin
        γ = 2.6752218744e8
        δg, Δ, σ, Gmax = 4.0e-3, 0.1, 0.9, 0.55
        D = 1.2                                     # ×10⁻¹⁰ m² s⁻¹
        k = (γ * δg * σ * Gmax)^2 * (Δ - δg / 3) * 1e-10
        gradients = collect(0.05:0.05:0.95)
        amps = [100 * exp(-k * g^2 * D) for g in gradients]

        ds = peakdataset(amps, :gradient, gradients; noise=0.2)
        expt = DiffusionExperiment(ds; γ, δ=δg, Δ, σ, Gmax, temp=298.15, solvent=:h2o,
                                   regions=signalregion())
        @test fitaxis(expt) == :gradient
        @test primaryparam(expt) == :D

        res = analyse1d(expt)
        @test Measurements.value(param(res[1], :D)) ≈ D rtol = 0.05

        η = viscosity(:h2o, 298.15)
        expected = 1.38e-23 * 298.15 / (6π * η * 0.001 * D * 1e-10) * 1e10
        @test param(res[1], :viscosity) ≈ η
        @test Measurements.value(param(res[1], :rH)) ≈ expected rtol = 0.05

        # with no solvent or temperature, D is still reported but rH is not derived
        dsn = peakdataset(amps, :gradient, gradients; noise=0.2)
        resn = analyse1d(DiffusionExperiment(dsn; γ, δ=δg, Δ, σ, Gmax,
                                             regions=signalregion()))
        @test !resn[1].postfitted
        @test_throws KeyError param(resn[1], :rH)

        @test shapefactor("SMSQ10.100") == 0.9
        @test shapefactor("SINE.100") ≈ 0.6366
        @test shapefactor("RECT") == 1.0
        @test shapefactor(nothing) == 1.0
        @test solventname("D2O") === :d2o
        @test solventname("H2O+D2O") === :h2o
        @test solventname("CDCl3") === nothing
    end

    @testset "Kinetics and series.csv" begin
        times = [0.0, 60.0, 120.0, 180.0]
        traces = [Trace(δ, lorentzian(100 - 0.3t, 2.0) .+ lorentzian(0.3t, 4.0))
                  for t in times]
        ds = Dataset1D(Planes(traces, [(; time=t) for t in times]), 8.0, "synthetic")
        expt = KineticsExperiment(ds;
                                  regions=[Region("reactant", 1.7, 2.3),
                                           Region("product", 3.7, 4.3)])
        res = analyse1d(expt)
        @test length(res) == 2
        @test isempty(res[1].parameters)                 # NoFitting fits nothing
        s1, s2 = only(res[1].series), only(res[2].series)
        @test s1.x == times
        @test length(s1.y) == 4
        # reactant falls, product rises
        @test Measurements.value(s1.y[1]) > Measurements.value(s1.y[end])
        @test Measurements.value(s2.y[1]) < Measurements.value(s2.y[end])

        # the measurements, long: one row per region per time point. This is the whole
        # deliverable of a kinetics run, which results.csv cannot carry because nothing
        # was fitted
        header, rows = seriestable(expt, ds, res)
        @test header == ["source", "label", "plane", "time (s)", "I", "I_err", "I_fit"]
        @test length(rows) == 8                          # 2 regions × 4 times
        @test all(length(row) == 7 for row in rows)
        @test all(row[1] == "synthetic" for row in rows)
        @test rows[1][2] == "reactant"
        @test [row[3] for row in rows[1:4]] == ["1", "2", "3", "4"]
        @test parse(Float64, rows[1][4]) == 0.0
        @test parse(Float64, rows[4][4]) == 180.0
        @test parse(Float64, rows[1][5]) > parse(Float64, rows[4][5])
        @test all(row[7] == "NA" for row in rows)        # NoFitting: no fitted curve

        rheader, rrows = resultstable(expt, res)
        @test rheader == ["label"]      # bounds are input, and live in regionlist.csv
        @test length(rrows) == 2
        @test occursin("Kinetics", experimentinfo(expt))
        @test occursin("Number of spectra: 4", experimentinfo(expt))
        @test occursin("NMRAnalysis.jl", experimentinfo(expt))
    end

    @testset "series.csv with several series" begin
        # a grouping variable becomes a key column in both files; long form means separate
        # runs need not share a fit axis
        traces = [peaktrace(100.0, 8.0) for _ in 1:4]
        vars = [(; time=0.0, run=1), (; time=1.0, run=1),
                (; time=0.0, run=2), (; time=2.0, run=2)]
        ds = Dataset1D(Planes(traces, vars), 0.0, "synthetic")
        expt = KineticsExperiment(ds; regions=signalregion())
        res = analyse1d(expt)
        @test length(res) == 1                  # one region, two runs
        @test length(only(res).series) == 2

        header, rows = seriestable(expt, ds, res)
        @test header ==
              ["source", "label", "plane", "run", "time (s)", "I", "I_err", "I_fit"]
        @test length(rows) == 4
        @test [row[4] for row in rows] == ["1", "1", "2", "2"]
        @test [parse(Float64, row[5]) for row in rows] == [0.0, 1.0, 0.0, 2.0]

        rheader, rrows = resultstable(expt, res)
        @test rheader == ["label"]              # NoFitting: nothing to report per region
        @test length(rrows) == 1
    end

    @testset "Column headers and one row per region" begin
        @test csvcolumn("R", "s-1") == "R (s-1)"
        @test csvcolumn("pB", "") == "pB"
        @test safename("L23N") == "L23N"
        @test safename("amide/1") == "amide_1"
        @test safename("") == "unnamed"

        times = [0.0, 0.005, 0.01, 0.02, 0.04]
        Rt, Ra = 16.0, 54.0
        f, ωN = tractf(; B0=18.8), 2π * 81.08e6
        traces = vcat([peaktrace(100 * exp(-Rt * t), 8.0) for t in times],
                      [peaktrace(100 * exp(-Ra * t), 8.0) for t in times])
        vars = vcat([(; time=t, which=:trosy) for t in times],
                    [(; time=t, which=:anti) for t in times])
        ds = Dataset1D(Planes(traces, vars), 0.0, "trosy.fid",
                       vcat(fill("trosy.fid", 5), fill("anti.fid", 5)))
        res = analyse1d(TractExperiment(ds; ωN, f, regions=signalregion()))
        expt = TractExperiment(ds; ωN, f, regions=signalregion())
        header, rows = resultstable(expt, res)
        @test header[1] == "label"
        @test "R_trosy (s-1)" in header && "R_trosy_err (s-1)" in header
        @test "R_anti (s-1)" in header && "R_anti_err (s-1)" in header
        @test "tauc (ns)" in header && "tauc_err (ns)" in header
        @test header[2] == "tauc (ns)"          # primaryparam comes first

        # one row per region, both components on it
        @test length(rows) == 1
        @test only(rows)[1] == "signal"
        @test all(cell != "NA" for cell in only(rows))

        # per-plane sources reach series.csv
        sheader, srows = seriestable(expt, ds, res)
        @test sheader ==
              ["source", "label", "plane", "which", "time (s)", "I", "I_err", "I_fit"]
        @test Set(row[1] for row in srows) == Set(["trosy.fid", "anti.fid"])
        @test all(row[end] != "NA" for row in srows)     # fitted curve at the data points
    end

    @testset "Result display" begin
        # a RegionResult prints as one line, not as a dump of every reduced point
        ds = peakdataset([100.0, 50.0, 25.0], :time, [0.0, 0.5, 1.0]; noise=0.5)
        res = analyse1d(RelaxationExperiment(ds; regions=signalregion()))
        text = sprint(show, res[1])
        @test count('\n', text) == 0
        @test startswith(text, "signal: ")
        @test occursin("R = ", text)
        @test length(sprint(show, res)) < 500
    end

    @testset "Region list" begin
        regs = [Region("signal", 7.9, 8.5), Region("wide", 1.0, 3.0)]
        header, rows = regionlisttable(regs, -1.0, 0.28)
        @test header == ["label", "lo (ppm)", "hi (ppm)"]
        @test length(rows) == 3                       # two regions plus the noise marker
        @test rows[1] == ["signal", "7.9", "8.5"]

        # the noise marker is written as a region as wide as the widest signal region, so
        # the row means something on its own; only its centre is read back
        noise = rows[end]
        @test noise[1] == NOISE_LABEL
        @test parse(Float64, noise[2]) ≈ -2.0
        @test parse(Float64, noise[3]) ≈ 0.0

        # with no regions at all it falls back to the default width
        _, narrow = regionlisttable(Region[], 5.0, 0.28)
        @test parse(Float64, narrow[end][3]) - parse(Float64, narrow[end][2]) ≈ 0.28

        # round trip
        dir = mktempdir()
        path = joinpath(dir, "regionlist.csv")
        open(path, "w") do f
            println(f, "# NMRAnalysis.jl test")
            println(f, join(header, ","))
            for row in rows
                println(f, join(row, ","))
            end
        end
        state = Dict{Symbol,Any}(:regions => Ref(Region[]), :active => Ref(0),
                                 :noisec => Ref(0.0))
        @test readregions!(state, path) == 2          # the noise row is not a region
        @test [r.label for r in state[:regions][]] == ["signal", "wide"]
        @test state[:regions][][1].lo ≈ 7.9
        @test state[:noisec][] ≈ -1.0

        # a results.csv from before the two files were separated still loads: label/lo/hi
        # columns and a noise comment
        legacy = joinpath(dir, "results.csv")
        open(legacy, "w") do f
            println(f, "# Noise position / ppm: 4.25")
            println(f, "region,lo,hi,A,A_err")
            println(f, "signal,7.9,8.5,1.0,0.1")
        end
        state2 = Dict{Symbol,Any}(:regions => Ref(Region[]), :active => Ref(0),
                                  :noisec => Ref(0.0))
        @test readregions!(state2, legacy) == 1
        @test state2[:noisec][] ≈ 4.25

        @test_throws ArgumentError readregions!(state2, joinpath(dir, "nope.csv"))
    end

    @testset "The reproduce line" begin
        call = analysiscall("relaxation1d", "11"; tau=[0.01, 0.03], model=:exponential)
        text = callstring(call, [Region("signal", 7.9, 8.5)], -1.0)
        @test occursin("relaxation1d(\"11\";", text)
        @test occursin("tau=[0.01, 0.03]", text)
        @test occursin("model=:exponential", text)
        # a single region is expressible as one integration triple, so the line repeats the
        # whole analysis without the window opening
        @test occursin("integration=(peakppm=8.2, noiseppm=-1.0, ppmwidth=0.6)", text)

        # several regions are not, so they are left out and the summary points at Load
        multi = callstring(call, [Region("a", 1.0, 2.0), Region("b", 3.0, 4.0)], -1.0)
        @test !occursin("integration", multi)

        # values with no literal are left out rather than printed misleadingly
        @test isempty(analysiscall("kinetics1d", "21"; model=NoFitting()).kwargs)
        @test callvalue(NoFitting()) === nothing
        @test callvalue(0.05:0.05:0.15) == repr([0.05, 0.1, 0.15])

        # a spectrum passed as an object has no literal either, so no line is printed at all
        @test isnothing(callstring(analysiscall("relaxation1d", Dict()), Region[], 0.0))

        # and it reaches summary.txt
        times = [0.0, 0.1, 0.2, 0.4]
        ds = peakdataset([100 * exp(-2.0t) for t in times], :time, times; noise=0.5)
        expt = RelaxationExperiment(ds; regions=signalregion())
        res = analyse1d(expt)
        path = joinpath(mktempdir(), "summary.txt")
        writesummary(path, expt, ds, res, expt.regions, call)
        written = read(path, String)
        @test occursin("Reproduce:", written)
        @test occursin("relaxation1d(", written)
        @test occursin("Relaxation rate", written)
        @test occursin("NMRAnalysis.jl", written)

        # without a recorded call there is simply no such section
        path2 = joinpath(mktempdir(), "summary.txt")
        writesummary(path2, expt, ds, res, expt.regions)
        @test !occursin("Reproduce:", read(path2, String))
    end

    @testset "Parameter prompts" begin
        # nothing is asked outside an interactive session: a default is taken, and a value
        # with no default raises an error naming what to supply
        @test ask("Diffusion delay Δ", 0.1; prompt=false) == 0.1
        @test_throws ArgumentError ask("Diffusion delay Δ"; prompt=false)
        @test_throws ArgumentError ask("Diffusion delay Δ", nothing; prompt=false)
        @test_throws ArgumentError askvector("relaxation delays", 4; prompt=false)
        @test askchoice("Which model?", ["a", "b"]; prompt=false) == 1
        @test askchoice("Which model?", ["a", "b"]; default=2, prompt=false) == 2

        @test parsevector("1 2 3") == [1.0, 2.0, 3.0]
        @test parsevector("1,2, 3") == [1.0, 2.0, 3.0]
        @test parsevector("1.5e-2 0.03") == [0.015, 0.03]
        @test parsevector("") === nothing
        @test parsevector("1 2 three") === nothing

        # a vdlist is a file of one value per line
        path = joinpath(mktempdir(), "vdlist")
        write(path, "0.01\n0.05\n0.1\n")
        @test parsevector(path) == [0.01, 0.05, 0.1]

        # acqus lookups answer `nothing` rather than throwing, so a precedence chain can
        # fall through them
        @test acqusvalue(nothing, :vdlist) === nothing
    end

    @testset "Results are written to disk" begin
        times = [0.0, 0.1, 0.2, 0.4]
        ds = peakdataset([100 * exp(-2.0t) for t in times], :time, times; noise=0.5)
        expt = RelaxationExperiment(ds; regions=signalregion())
        res = analyse1d(expt)

        dir = mktempdir()
        writeresults!(expt, ds, res, expt.regions, dir)

        @test isfile(joinpath(dir, "results.csv"))
        @test isfile(joinpath(dir, "series.csv"))
        @test isfile(joinpath(dir, "regionlist.csv"))
        @test isfile(joinpath(dir, "regions", "signal.csv"))

        datalines(f) = filter(l -> !startswith(l, "#"), readlines(joinpath(dir, f)))
        lines = datalines("series.csv")
        @test length(lines) == 5                       # header plus one row per spectrum
        @test split(lines[1], ',') ==
              ["source", "label", "plane", "time (s)", "I", "I_err", "I_fit"]
        @test datalines(joinpath("regions", "signal.csv")) == lines
        @test occursin("Noise position", read(joinpath(dir, "results.csv"), String))

        # every quantity a diffusion analysis reports describes its region, the solvent
        # viscosity included, so they all sit on the one row
        γ, δg, Δ, σ, Gmax = 2.6752218744e8, 4.0e-3, 0.1, 0.9, 0.55
        k = (γ * δg * σ * Gmax)^2 * (Δ - δg / 3) * 1e-10
        g = collect(0.05:0.05:0.95)
        dsd = peakdataset([100 * exp(-k * gi^2 * 1.2) for gi in g], :gradient, g;
                          noise=0.2)
        exptd = DiffusionExperiment(dsd; γ, δ=δg, Δ, σ, Gmax, temp=298.15, solvent=:h2o,
                                    regions=signalregion())
        resd = analyse1d(exptd)

        dird = mktempdir()
        writeresults!(exptd, dsd, resd, exptd.regions, dird)

        rheader, rrows = resultstable(exptd, resd)
        @test "D (1e-10 m2/s)" in rheader
        @test "rH (A)" in rheader
        @test "viscosity (mPa s)" in rheader
        @test length(rrows) == 1
    end
end
