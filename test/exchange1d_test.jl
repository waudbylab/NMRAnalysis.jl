using NMRAnalysis
import NMRAnalysis.Exchange1D:
    NoExchangeModel, TwoStateModel, TwoStateBindingModel,
    R1Experiment, CESTExperiment, ExchangeProblem,
    nstates, modelname, nmolecules,
    exchange_matrix, populations, default_params,
    default_spin_params, default_nuisance_params,
    field_label, liouvillian,
    simulate!, residuals, fit
using ComponentArrays
using Measurements
using LinearAlgebra
using Test

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

        params = ComponentArray(; model=default_params(model))
        K = exchange_matrix(model, params, Dict{String,Float64}())
        @test size(K) == (1, 1)
        @test K[1, 1] == 0.0

        p = populations(model, params, Dict{String,Float64}())
        @test p == [1.0]
    end

    @testset "TwoStateModel" begin
        model = TwoStateModel()
        @test nstates(model) == 2

        params = ComponentArray(;
            model=ComponentArray(; kex=1000.0, pB=0.05),
        )
        K = exchange_matrix(model, params, Dict{String,Float64}())
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
        model = TwoStateBindingModel(Dict(:X => "protein", :Y => "ligand"))
        @test nstates(model) == 2
        @test nmolecules(model) == 2

        params = ComponentArray(;
            model=ComponentArray(; Kd=100.0, koff=5000.0),
        )
        conc = Dict("protein" => 100.0, "ligand" => 200.0)

        K = exchange_matrix(model, params, conc)
        @test size(K) == (2, 2)
        @test sum(K; dims=1) ≈ zeros(1, 2) atol = 1e-10

        p = populations(model, params, conc)
        @test sum(p) ≈ 1.0
        @test 0 < p[2] < 1  # some fraction bound

        # empty moleculemap should error
        model_empty = TwoStateBindingModel()
        @test_throws ArgumentError exchange_matrix(model_empty, params, conc)
    end

    @testset "Binding fraction" begin
        # when Kd = 0, all should be bound (pB = min(Yt, Xt) / Xt)
        # when Kd >> concentrations, pB ≈ 0

        # equal concentrations, Kd = 0 limit (use very small Kd)
        pB = Exchange1D._binding_fraction(1e-10, 100.0, 100.0)
        @test pB ≈ 1.0 atol = 1e-4

        # Kd much larger than concentrations
        pB = Exchange1D._binding_fraction(1e6, 100.0, 100.0)
        @test pB < 0.01

        # symmetric case: Kd = Xt = Yt
        pB = Exchange1D._binding_fraction(100.0, 100.0, 100.0)
        # quadratic gives pB ≈ 0.382
        @test 0.3 < pB < 0.5
    end

    @testset "Liouvillian - NoExchange" begin
        model = NoExchangeModel()
        params = ComponentArray(;
            model=ComponentArray(),
            spin=ComponentArray(;
                delta=[0.0],
                R2_14p1T=[15.0],
                R1_14p1T=[1.5],
            ),
        )

        L = liouvillian(model, params, 14.1, 600e6, 0.0, 500.0, Dict{String,Float64}())
        @test size(L) == (3, 3)

        # diagonal: -R2, -R2, -R1
        @test L[1, 1] ≈ -15.0
        @test L[2, 2] ≈ -15.0
        @test L[3, 3] ≈ -1.5

        # spin-lock couples My ↔ Mz
        @test L[2, 3] ≈ -2π * 500.0
        @test L[3, 2] ≈ 2π * 500.0
    end

    @testset "Liouvillian - TwoState" begin
        model = TwoStateModel()
        N = nstates(model)
        params = ComponentArray(;
            model=ComponentArray(; kex=1000.0, pB=0.05),
            spin=ComponentArray(;
                delta=[-62.0, -58.0],
                R2_14p1T=[15.0, 150.0],
                R1_14p1T=[1.5],
            ),
        )

        L = liouvillian(model, params, 14.1, 600e6, -60.0 * 600.0, 500.0,
                        Dict{String,Float64}())
        @test size(L) == (3N, 3N)
        @test size(L) == (6, 6)

        # R1 should be the same for both states (shared, length-1 vector)
        @test L[3, 3] ≈ L[6, 6]  # both Mz diagonal entries have same R1

        # R2 should differ between states
        @test L[1, 1] != L[4, 4]  # Mx diagonals differ

        # exchange terms should be present in off-diagonal blocks
        # K[2,1] = kex * pB = 50, added to L[4,1], L[5,2], L[6,3]
        @test L[4, 1] ≈ 1000.0 * 0.05   # Mx exchange A→B
        @test L[5, 2] ≈ 1000.0 * 0.05   # My exchange A→B
        @test L[6, 3] ≈ 1000.0 * 0.05   # Mz exchange A→B
    end

    @testset "Liouvillian eigenvalues" begin
        # for no exchange, eigenvalues should be related to R1, R2
        model = NoExchangeModel()
        params = ComponentArray(;
            model=ComponentArray(),
            spin=ComponentArray(;
                delta=[0.0],
                R2_14p1T=[15.0],
                R1_14p1T=[1.5],
            ),
        )

        # on-resonance with no spin-lock: eigenvalues = -R2, -R2, -R1
        L = liouvillian(model, params, 14.1, 600e6, 0.0, 0.0, Dict{String,Float64}())
        eigenvals = sort(real.(eigvals(L)))
        @test eigenvals[1] ≈ -15.0 atol = 1e-10
        @test eigenvals[2] ≈ -15.0 atol = 1e-10
        @test eigenvals[3] ≈ -1.5 atol = 1e-10
    end

    @testset "R1 simulate! - exponential decay" begin
        delays = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
        observed = [1.0 ± 0.02, 0.9 ± 0.02, 0.8 ± 0.02,
                    0.6 ± 0.02, 0.4 ± 0.02, 0.2 ± 0.02]
        predicted = zeros(length(delays))
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(), delays,
                            observed, predicted, :exponential_decay)

        model = NoExchangeModel()
        params = ComponentArray(;
            model=ComponentArray(),
            spin=ComponentArray(;
                R2_14p1T=[15.0],
                R1_14p1T=[2.0],
            ),
            nuisance=ComponentArray(;
                R1_14p1T_I0=1.0,
            ),
        )

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
        params = ComponentArray(;
            model=ComponentArray(),
            spin=ComponentArray(;
                R2_14p1T=[15.0],
                R1_14p1T=[2.0],
            ),
            nuisance=ComponentArray(;
                R1_14p1T_I0=1.0,
                R1_14p1T_inv_factor=2.0,
            ),
        )

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

    @testset "default_params - structure" begin
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                            [0.1, 0.2, 0.5],
                            [1.0 ± 0.02, 0.8 ± 0.02, 0.5 ± 0.02],
                            zeros(3), :exponential_decay)

        prob = ExchangeProblem([expt], TwoStateModel())
        params = default_params(prob)

        # model section
        @test haskey(params, :model)
        @test params.model.kex == 1000.0
        @test params.model.pB == 0.05

        # spin section — R1 experiments don't need delta (chemical shifts)
        @test haskey(params, :spin)
        @test !haskey(params.spin, :delta)

        fl = field_label(14.1)
        @test haskey(params.spin, Symbol("R2_", fl))
        @test haskey(params.spin, Symbol("R1_", fl))
        @test length(params.spin[Symbol("R2_", fl)]) == 2  # per state
        @test length(params.spin[Symbol("R1_", fl)]) == 1  # scalar (shared across states)

        # nuisance section — flat keys like :R1_14p1T_I0
        @test haskey(params, :nuisance)
        @test haskey(params.nuisance, Symbol("R1_", fl, "_I0"))
        @test params.nuisance[Symbol("R1_", fl, "_I0")] == 1.0
    end

    @testset "default_params - inversion recovery" begin
        expt = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                            [0.1, 0.2, 0.5],
                            [1.0 ± 0.02, 0.8 ± 0.02, 0.5 ± 0.02],
                            zeros(3), :inversion_recovery)

        prob = ExchangeProblem([expt], NoExchangeModel())
        params = default_params(prob)

        fl = field_label(14.1)
        @test params.nuisance[Symbol("R1_", fl, "_I0")] == 1.0
        @test params.nuisance[Symbol("R1_", fl, "_inv_factor")] == 2.0
    end

    @testset "default_params - multiple fields" begin
        expt1 = R1Experiment(nothing, 14.1, Dict{String,Float64}(),
                             [0.1], [1.0 ± 0.02], zeros(1),
                             :exponential_decay)
        expt2 = R1Experiment(nothing, 18.79, Dict{String,Float64}(),
                             [0.1], [1.0 ± 0.02], zeros(1),
                             :exponential_decay)

        prob = ExchangeProblem([expt1, expt2], TwoStateModel())
        params = default_params(prob)

        fl1 = field_label(14.1)
        fl2 = field_label(18.79)
        @test haskey(params.spin, Symbol("R1_", fl1))
        @test haskey(params.spin, Symbol("R1_", fl2))
        @test haskey(params.spin, Symbol("R2_", fl1))
        @test haskey(params.spin, Symbol("R2_", fl2))

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

        params = ComponentArray(;
            model=ComponentArray(),
            spin=ComponentArray(;
                R2_14p1T=[15.0],
                R1_14p1T=[R1_true],
            ),
            nuisance=ComponentArray(;
                R1_14p1T_I0=I0_true,
            ),
        )

        # simulate should fill predicted_intensities
        simulate!(prob, params)
        @test expt.predicted_intensities ≈ true_intensities

        # residuals should be near zero when predicted matches observed
        r = residuals(prob, params)
        @test all(abs.(r) .< 1e-10)
    end
end
