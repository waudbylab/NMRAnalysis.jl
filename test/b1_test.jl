# B₁ calibration and inhomogeneity (src/b1.jl): the power-to-field law, its inverse, and
# the quadrature the simulations average over. No data files and no spectra: a calibration
# is a handful of numbers.

using NMRAnalysis
using NMRAnalysis: refpower, ν1ref, npoints, source, calibrationfit, B1NODES, B1WEIGHTS,
                   DEFAULTINHOMOGENEITY
using NMRTools: Power, db, watts, hz, nucleus
using Measurements
using Test

@testset "B1Calibration" begin
    @testset "one point is the ideal power law" begin
        # The nominal calibration a spectrum's own reference pulse gives, built from the
        # same numbers NMRTools' four-argument `hz` takes. The two must agree exactly, or
        # simply adopting the calibration would shift every spin-lock field in an
        # already-published analysis.
        pulse, refp = 13.29e-6, Power(-9.03, :dB)
        cal = B1Calibration([refp], [90 / (360 * pulse)])

        @test npoints(cal) == 1
        @test Measurements.value(linearity(cal)) == 1.0
        @test inhomogeneity(cal) == DEFAULTINHOMOGENEITY
        for dB in (-9.03, 0.0, 6.0, 20.0)
            p = Power(dB, :dB)
            @test hz(p, cal) === hz(p, refp, pulse, 90)
        end
        # at the reference power the field is the measured one, not one round-tripped
        # through a logarithm
        @test hz(refp, cal) === 90 / (360 * pulse)
        # 6 dB more attenuation roughly halves the field (exactly, a factor 10^-0.3)
        @test hz(Power(db(refp) + 6, :dB), cal) ≈ hz(refp, cal) * 10^(-0.3)
        @test hz(Power(db(refp) + 6, :dB), cal) ≈ hz(refp, cal) / 2 rtol = 5e-3
    end

    @testset "several points fit the amplifier's response" begin
        # A response 2% short of the ideal law, sampled 6 dB apart
        dB = [-18.0, -12.0, -6.0]
        L, ν0 = 0.98, 500.0
        ν = [ν0 * 10^(-L * (p - dB[1]) / 20) for p in dB]
        cal = B1Calibration(Power.(dB, :dB), ν .± 1.0; nuc=nucleus("19F"),
                            inhomogeneity=0.043, source="synthetic")

        @test npoints(cal) == 3
        @test db(refpower(cal)) == dB[1]
        @test Measurements.value(ν1ref(cal)) ≈ ν0 rtol = 1e-6
        @test Measurements.value(linearity(cal)) ≈ L rtol = 1e-6
        @test Measurements.uncertainty(linearity(cal)) > 0
        @test inhomogeneity(cal) == 0.043
        @test nucleus(cal) == nucleus("19F")
        @test occursin("synthetic", source(cal))

        # the curve passes through the measurements, and interpolates between them
        for i in eachindex(dB)
            @test hz(Power(dB[i], :dB), cal) ≈ ν[i] rtol = 1e-6
        end
        @test ν[2] < hz(Power(-15.0, :dB), cal) < ν[1]   # dB is attenuation

        # `Power` inverts `hz`: the power to set for a wanted spin-lock strength
        for target in (250.0, 1000.0, 4000.0)
            @test hz(Power(target, cal), cal) ≈ target rtol = 1e-9
        end
        @test db(Power(ν[2], cal)) ≈ dB[2] rtol = 1e-9
        # and a shortfall means more power than the ideal law would ask for
        @test watts(Power(2ν0, cal)) > watts(Power(dB[1] - 6, :dB))
    end

    @testset "a calibration is one object, not a list" begin
        cal = B1Calibration([Power(-12.0, :dB)], [1000.0])
        @test hz.([Power(-12.0, :dB), Power(-6.0, :dB)], cal) ≈ [1000.0, 501.187] rtol = 1e-4
    end

    @testset "what cannot be calibrated" begin
        @test_throws ArgumentError B1Calibration(Power{Float64}[], Float64[])
        @test_throws ArgumentError B1Calibration([Power(-12.0, :dB)], [1000.0, 2000.0])
        @test_throws ArgumentError B1Calibration([Power(-12.0, :dB)], [1000.0];
                                                 inhomogeneity=5)
        @test_throws ArgumentError B1Calibration([Power(-12.0, :dB)], [-1000.0])
        # two measurements at the same power say nothing about the slope
        @test_throws ArgumentError B1Calibration(Power.([-12.0, -12.0], :dB),
                                                 [1000.0, 1010.0])
        @test_throws ArgumentError Power(-5.0, B1Calibration([Power(-12.0, :dB)], [1000.0]))
    end

    @testset "the fit is weighted by the uncertainties" begin
        # One point pulled off the line, first with a large uncertainty and then with a
        # small one: the fitted slope should follow the precise measurement.
        dB = [-18.0, -12.0, -6.0]
        ν = [500.0, 1000.0, 2000.0]
        loose = calibrationfit(dB, [ν[1] ± 1.0, (ν[2] * 1.2) ± 500.0, ν[3] ± 1.0], dB[1])
        tight = calibrationfit(dB, [ν[1] ± 1.0, (ν[2] * 1.2) ± 1.0, ν[3] ± 1.0], dB[1])
        @test Measurements.value(loose[1]) ≈ 500.0 rtol = 0.02
        @test Measurements.uncertainty(tight[2]) < Measurements.uncertainty(loose[2])
    end
end

@testset "B1Distribution" begin
    @testset "moments of the three-point quadrature" begin
        σ = 0.05
        d = B1Distribution(σ)
        # Gauss-Hermite: exact for the mean, variance and fourth moment of a Gaussian, for
        # any σ. This is what fixed node positions cannot do - which is the whole reason
        # the nodes scale with σ.
        @test sum(d.weight) ≈ 1.0
        @test sum(d.weight .* d.scaling) ≈ 1.0
        @test sum(d.weight .* (d.scaling .- 1) .^ 2) ≈ σ^2
        @test sum(d.weight .* (d.scaling .- 1) .^ 4) ≈ 3σ^4
        @test d.scaling ≈ [1 - √3 * σ, 1.0, 1 + √3 * σ]
        @test collect(d.weight) ≈ collect(B1WEIGHTS)
        @test length(B1NODES) == 3
    end

    @testset "no inhomogeneity is one field" begin
        d = B1Distribution(0.0)
        @test d.scaling == [1.0]
        @test d.weight == [1.0]
        @test b1average(ν -> 2ν, d, 500.0) == 1000.0
    end

    @testset "from a calibration, and by hand" begin
        cal = B1Calibration([Power(-12.0, :dB)], [1000.0]; inhomogeneity=0.08)
        @test B1Distribution(cal).scaling ≈ B1Distribution(0.08).scaling
        # any other distribution, or a profile measured on the probe, is just nodes and
        # weights
        d = B1Distribution([0.9, 1.0, 1.2], [1.0, 2.0, 1.0])
        @test sum(d.weight) ≈ 1.0
        @test d.weight ≈ [0.25, 0.5, 0.25]
        @test_throws ArgumentError B1Distribution([0.0, 1.0], [1.0, 1.0])
        @test_throws ArgumentError B1Distribution([0.9, 1.1], [1.0])
        @test_throws ArgumentError B1Distribution(-0.01)
        @test_throws ArgumentError B1Distribution(0.6)      # 1 - √3σ ≤ 0
    end

    @testset "averaging an observable" begin
        d = B1Distribution(0.05)
        ν̄ = 500.0
        # linear in the field: the average is the nominal field
        @test b1average(ν -> ν, d, ν̄) ≈ ν̄
        # quadratic: ⟨ν²⟩ = ν̄²(1 + σ²), which the quadrature gets exactly
        @test b1average(ν -> ν^2, d, ν̄) ≈ ν̄^2 * (1 + 0.05^2)
    end

    @testset "the decay a spread of fields gives" begin
        d = B1Distribution(0.05)
        rate(ν) = 30 / (1 + (ν / 300)^2)
        ν̄ = 400.0
        t = [0.0, 0.02, 0.05, 0.1]

        # the weighted sum of exponentials, written out
        expected = [sum(w * exp(-rate(s * ν̄) * tᵢ)
                        for (s, w) in zip(d.scaling, d.weight)) for tᵢ in t]
        @test b1decay(rate, d, ν̄, t) ≈ expected
        # normalised at t = 0 whatever the spread
        @test b1decay(rate, d, ν̄, [0.0]) ≈ [1.0]
        # with no inhomogeneity it is the single exponential it was before
        @test b1decay(rate, B1Distribution(0.0), ν̄, t) ≈ exp.(-rate(ν̄) .* t)

        # a sum of exponentials decays more slowly than the exponential of the mean rate:
        # the slow-relaxing part of the sample dominates what is left at long times
        mono = exp.(-b1average(rate, d, ν̄) .* t)
        @test all(b1decay(rate, d, ν̄, t)[2:end] .> mono[2:end])

        # b1rate is the rate of the single exponential through that decay at whatever
        # time it is matched at
        T = 0.1
        @test only(b1decay(rate, d, ν̄, (T,))) ≈ exp(-b1rate(rate, d, ν̄, T) * T)
    end

    @testset "averaging a fitted rate" begin
        d = B1Distribution(0.05)
        # A rate that varies steeply with the field, as R₁ρ does at low spin-lock
        rate(ν) = 30 / (1 + (ν / 300)^2)
        ν̄, T = 400.0, 0.1
        mean_rate = b1average(rate, d, ν̄)

        # what a monoexponential fit to the averaged decay reports is lower than the mean
        # of the rates (the decay is a sum of exponentials, and Jensen's inequality bites)
        @test b1rate(rate, d, ν̄, T) < mean_rate
        # and approaches it as the evolution time shrinks
        @test b1rate(rate, d, ν̄, 1e-6) ≈ mean_rate rtol = 1e-4
        @test b1rate(rate, d, ν̄, 0.0) == mean_rate
        # with no inhomogeneity there is nothing to average
        @test b1rate(rate, B1Distribution(0.0), ν̄, T) ≈ rate(ν̄)
    end
end
