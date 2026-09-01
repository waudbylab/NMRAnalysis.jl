using NMRAnalysis
import NMRAnalysis.Exchange1D:
                               NoExchangeModel, TwoStateModel, TwoStateBindingModel,
                               ThreeStateModel, ThreeStateBindingModel,
                               R1Experiment, CESTExperiment, ExchangeProblem,
                               nstates, modelname, nmolecules,
                               exchangematrix, populations, defaultparams,
                               default_spin_params, default_nuisance_params,
                               field_label, liouvillian, liouvillian_inhom,
                               AbstractExperiment,
                               simulate!, residuals, fit,
                               sampleconcentrations, moleculeconcentration, experimentgroups
using ComponentArrays
using Measurements
using Test

# Minimal experiment stand-in for tests that build a Liouvillian directly:
# liouvillian/liouvillian_inhom only reach for `field_teslas` (via field_label)
# and the base frequency, `spec[1, :bf]`.
struct StubSpec
    bf::Float64
end
Base.getindex(s::StubSpec, ::Int, ::Symbol) = s.bf

struct StubExperiment <: AbstractExperiment
    field_teslas::Float64
    spec::StubSpec
end

# Minimal spec stand-in for tests that need `short_expt_path` (via `expt.spec[:filename]`),
# e.g. to check the wording of an error message.
struct StubPathSpec
    filename::String
end
Base.getindex(s::StubPathSpec, ::Symbol) = s.filename

@testset "Exchange1D" begin
    @testset "field_label" begin
        @test field_label(14.1) == Symbol("14p1T")
        @test field_label(23.49) == Symbol("23p49T")
        @test field_label(11.0) == Symbol("11p0T")
    end

    @testset "NoExchangeModel" begin
        model = NoExchangeModel()
        @test nstates(model) == 1
        @test modelname(model) == "No exchange"

        params = ComponentArray(; model=defaultparams(model))
        K = exchangematrix(model, params, Dict{String,Float64}())
        @test size(K) == (1, 1)
        @test K[1, 1] == 0.0

        p = populations(model, params, Dict{String,Float64}())
        @test p == [1.0]
    end

    @testset "TwoStateModel" begin
        model = TwoStateModel()
        @test nstates(model) == 2

        pB = 0.05
        dGB = log((1 - pB) / pB)  # ΔG reconstructing pB = 0.05
        params = ComponentArray(;
                                model=ComponentArray(; logkex=log(1000.0), dGB=dGB),)
        K = exchangematrix(model, params, Dict{String,Float64}())
        @test size(K) == (2, 2)

        # column sums must be zero (conservation of magnetisation)
        @test sum(K; dims=1) ≈ zeros(1, 2) atol = 1e-10

        # off-diagonal rates
        @test K[2, 1] ≈ 1000.0 * 0.05   # kBA = kex * pB
        @test K[1, 2] ≈ 1000.0 * 0.95   # kAB = kex * pA

        p = populations(model, params, Dict{String,Float64}())
        @test sum(p) ≈ 1.0
        @test p[1] ≈ 0.95
        @test p[2] ≈ 0.05
    end

    @testset "TwoStateBindingModel" begin
        model = TwoStateBindingModel(Dict(:A => "protein", :X => "ligand"))
        @test nstates(model) == 2
        @test nmolecules(model) == 2

        @test defaultparams(model).logKd ≈ log(100.0)
        @test defaultparams(model).logkoff ≈ log(5000.0)
    end

    @testset "modelname - descriptive parameter lists (#37)" begin
        @test modelname(NoExchangeModel()) == "No exchange"
        @test modelname(TwoStateModel()) == "2-state exchange (kex, pB)"
        @test modelname(TwoStateBindingModel()) == "2-state binding (Kd, koff)"
        @test modelname(ThreeStateModel()) == "3-state exchange (koffB, pB, koffC, pC)"
        @test modelname(ThreeStateBindingModel()) == "3-state binding (Kd1, koff1, Kd2, koff2)"
    end

    @testset "sampleconcentrations reads the stored field" begin
        # previously this recomputed from `expt.spec`, so an experiment built
        # without a real NMRData spec (as in these unit tests) would error
        expt = R1Experiment(nothing, 14.1, Dict("protein" => 50.0), [0.1],
                            [1.0 ± 0.02], zeros(1), :exponential_decay)
        @test sampleconcentrations(expt) == Dict("protein" => 50.0)
    end

    @testset "moleculeconcentration - clear error for unmatched sample (#36)" begin
        model = TwoStateBindingModel(Dict(:A => "protein", :X => "ligand"))
        spec = StubPathSpec("/data/sophia_trypsin/34/pdata/1")

        # neither sample metadata nor a fallback concentration is available
        unmatched = R1Experiment(spec, 14.1, Dict{String,Float64}(), [0.1],
                                 [1.0 ± 0.02], zeros(1), :exponential_decay)
        @test_throws ArgumentError moleculeconcentration(model, unmatched, :A)
        try
            moleculeconcentration(model, unmatched, :A)
        catch err
            @test err isa ArgumentError
            @test occursin("protein", err.msg)
            @test occursin(":A", err.msg)
            @test occursin("sophia_trypsin", err.msg)
        end

        # a fallback concentration (as entered interactively) resolves it
        model.concentrations["protein"] = 100.0
        @test moleculeconcentration(model, unmatched, :A) == 100.0

        # sample metadata, where present, takes precedence over the fallback
        matched = R1Experiment(spec, 14.1, Dict("protein" => 250.0), [0.1],
                               [1.0 ± 0.02], zeros(1), :exponential_decay)
        @test moleculeconcentration(model, matched, :A) == 250.0
    end

    @testset "experimentgroups - preserves first-seen order (#39)" begin
        groups = experimentgroups(iseven, [1, 3, 2, 4, 5])
        @test length(groups) == 2
        @test groups[1] == [1, 3, 5]  # odd, first seen
        @test groups[2] == [2, 4]     # even, first seen second

        singletons = experimentgroups(identity, [1, 2, 3])
        @test length(singletons) == 3
    end

    # Note: most liouvillian tests require a real NMR experiment (spec with :bf metadata)
    # and cannot be tested with spec=nothing. These are tested via integration tests.
    # The stub below supplies just enough (field_teslas, spec[1, :bf]) to build a
    # Liouvillian directly.

    @testset "liouvillian_inhom - equilibrium is stationary" begin
        # The exchange block must be K ⊗ I₃, i.e. K[i,j] = rate from j to i, matching
        # the convention documented on exchangematrix (zero column sums). Building it
        # transposed leaves the diagonal untouched, so it is easy to miss: it shows up
        # as magnetisation that is not conserved.
        #
        # Physical invariant: with the saturation field off, magnetisation starting at
        # equilibrium (Mz = populations) must not evolve, so L * M₀ = 0.
        model = TwoStateModel()
        pB = 0.05
        params = ComponentArray(;
                                model=ComponentArray(; logkex=log(1000.0),
                                                     dGB=log((1 - pB) / pB)),
                                spin=ComponentArray(; delta=[0.0, 5.0],
                                                    R1_14p1T=[1.0], R2_14p1T=[10.0, 10.0]))
        expt = StubExperiment(14.1, StubSpec(564.0))

        p0 = populations(model, params, expt)
        N = nstates(model)
        M0 = zeros(3N + 1)
        for i in 1:N
            M0[3(i - 1) + 3] = p0[i]
        end
        M0[end] = 1.0

        # spinlock_hz = 0: no saturation, so equilibrium must be a fixed point
        L = liouvillian_inhom(model, params, expt, -100.0, 0.0)
        @test L * M0 ≈ zeros(3N + 1) atol = 1e-9

        # the exchange block itself must conserve magnetisation (zero column sums)
        # over the Mz rows, independent of relaxation/shift terms
        K = exchangematrix(model, params, expt)
        @test sum(K; dims=1) ≈ zeros(1, N) atol = 1e-10
        @test K * p0 ≈ zeros(N) atol = 1e-10

        # and the homogeneous Liouvillian must agree with it on the exchange block
        Lhom = liouvillian(model, params, expt, -100.0, 0.0)
        @test Lhom[1:(3N), 1:(3N)] ≈ L[1:(3N), 1:(3N)] atol = 1e-10
    end

    @testset "R1 simulate! - exponential decay" begin
        delays = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
        observed = [1.0 ± 0.02, 0.9 ± 0.02, 0.8 ± 0.02,
                    0.6 ± 0.02, 0.4 ± 0.02, 0.2 ± 0.02]
        predicted = zeros(length(delays))
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(), delays,
                            observed, predicted, :exponential_decay)

        model = NoExchangeModel()
        params = ComponentArray(; model=ComponentArray(),
                                spin=ComponentArray(; R2_14p1T=[15.0],
                                                    R1_14p1T=[2.0],),
                                nuisance=ComponentArray(;
                                                        R1_14p1T_I0=1.0,),)

        simulate!(expt, model, params)

        R1 = 2.0
        expected = 1.0 .* exp.(-delays .* R1)
        @test expt.predicted_intensities ≈ expected
    end

    @testset "R1 simulate! - inversion recovery" begin
        delays = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
        observed = [-0.8 ± 0.02, -0.5 ± 0.02, 0.0 ± 0.02,
                    0.5 ± 0.02, 0.8 ± 0.02, 0.95 ± 0.02]
        predicted = zeros(length(delays))
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(), delays,
                            observed, predicted, :inversion_recovery)

        model = NoExchangeModel()
        params = ComponentArray(; model=ComponentArray(),
                                spin=ComponentArray(; R2_14p1T=[15.0],
                                                    R1_14p1T=[2.0],),
                                nuisance=ComponentArray(; R1_14p1T_I0=1.0,
                                                        R1_14p1T_inv_factor=2.0,),)

        simulate!(expt, model, params)

        R1 = 2.0
        expected = 1.0 .* (1.0 .- 2.0 .* exp.(-delays .* R1))
        @test expt.predicted_intensities ≈ expected
    end

    @testset "Residuals" begin
        delays = [0.1, 0.2, 0.5]
        noise = 0.05
        observed = [0.8 ± noise, 0.6 ± noise, 0.3 ± noise]
        predicted = [0.82, 0.61, 0.28]
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(), delays,
                            observed, predicted, :exponential_decay)

        r = residuals(expt)
        @test length(r) == 3
        @test r[1] ≈ (0.8 - 0.82) / noise
        @test r[2] ≈ (0.6 - 0.61) / noise
        @test r[3] ≈ (0.3 - 0.28) / noise
    end

    @testset "defaultparams - structure" begin
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                            [0.1, 0.2, 0.5],
                            [1.0 ± 0.02, 0.8 ± 0.02, 0.5 ± 0.02],
                            zeros(3), :exponential_decay)

        prob = ExchangeProblem([expt], TwoStateModel())
        params = defaultparams(prob)

        # model section
        @test haskey(params, :model)
        @test params.model.logkex ≈ log(1000.0)
        @test params.model.dGB ≈ log(0.95 / 0.05)

        # spin section — R1 experiments don't need delta (chemical shifts)
        @test haskey(params, :spin)
        @test !haskey(params.spin, :delta)

        fl = field_label(14.1)
        @test haskey(params.spin, Symbol("R1_", fl))
        @test length(params.spin[Symbol("R1_", fl)]) == 1  # scalar (shared across states)

        # nuisance section — flat keys like :R1_14p1T_I0
        @test haskey(params, :nuisance)
        @test haskey(params.nuisance, Symbol("R1_", fl, "_I0"))
        @test params.nuisance[Symbol("R1_", fl, "_I0")] == 1.0
    end

    @testset "defaultparams - inversion recovery" begin
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                            [0.1, 0.2, 0.5],
                            [1.0 ± 0.02, 0.8 ± 0.02, 0.5 ± 0.02],
                            zeros(3), :inversion_recovery)

        prob = ExchangeProblem([expt], NoExchangeModel())
        params = defaultparams(prob)

        fl = field_label(14.1)
        @test params.nuisance[Symbol("R1_", fl, "_I0")] == 1.0
        @test params.nuisance[Symbol("R1_", fl, "_inv_factor")] == 2.0
    end

    @testset "defaultparams - multiple fields" begin
        expt1 = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                             [0.1], [1.0 ± 0.02], zeros(1),
                             :exponential_decay)
        expt2 = R1Experiment(nothing, 18.79, Dict{String,Float64}(),
                             [0.1], [1.0 ± 0.02], zeros(1),
                             :exponential_decay)

        prob = ExchangeProblem([expt1, expt2], TwoStateModel())
        params = defaultparams(prob)

        fl1 = field_label(14.1)
        fl2 = field_label(18.79)
        @test haskey(params.spin, Symbol("R1_", fl1))
        @test haskey(params.spin, Symbol("R1_", fl2))

        # nuisance params for both fields
        @test haskey(params.nuisance, Symbol("R1_", fl1, "_I0"))
        @test haskey(params.nuisance, Symbol("R1_", fl2, "_I0"))
    end

    @testset "Problem simulate! and residuals" begin
        delays = [0.1, 0.2, 0.5, 1.0]
        R1_true = 2.0
        I0_true = 1.0
        true_intensities = I0_true .* exp.(-delays .* R1_true)
        noise = 0.02

        observed = true_intensities .± noise
        predicted = zeros(length(delays))
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(), delays,
                            observed, predicted, :exponential_decay)

        prob = ExchangeProblem([expt], NoExchangeModel())

        params = ComponentArray(; model=ComponentArray(),
                                spin=ComponentArray(; R2_14p1T=[15.0],
                                                    R1_14p1T=[R1_true],),
                                nuisance=ComponentArray(;
                                                        R1_14p1T_I0=I0_true,),)

        # simulate should fill predicted_intensities
        simulate!(prob, params)
        @test expt.predicted_intensities ≈ true_intensities

        # residuals should be near zero when predicted matches observed
        r = residuals(prob, params)
        @test all(abs.(r) .< 1e-10)
    end
end
