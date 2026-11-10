# GUI2D: the entry points' arguments, with no data and no window.

using NMRAnalysis
using NMRAnalysis.GUI2D: planelist, asexptpath, asksaturation
using NMRAnalysis: analysiscall, callstring, literal
using Test

@testset "GUI2D" begin
    @testset "Spectra given as paths or experiment numbers" begin
        @test planelist("1/pdata/1") == ["1/pdata/1"]
        @test planelist(3) == ["3"]
        @test planelist(10:12) == ["10", "11", "12"]
        @test planelist(["a", "b"]) == ["a", "b"]
        @test asexptpath(1:2) == ["1", "2"]
    end

    @testset "Data are positional, roles in order" begin
        @test hasmethod(hetnoe2d, Tuple{Any})                 # planes, with saturation
        @test hasmethod(hetnoe2d, Tuple{Any,Any})             # reference, saturated
        @test hasmethod(rdc2d, Tuple{Any,Any})                # isotropic, aligned
        @test hasmethod(ccr2d, Tuple{Any,Any})                # decay, buildup
        @test hasmethod(methylccr2d, Tuple{Any,Any})          # decay, buildup
        @test !hasmethod(relaxation2d, Tuple{Any,Any})        # delays are a keyword
        @test_throws ArgumentError hetnoe2d(["a"], ["b", "c"])
    end

    @testset "Times are T, and former names still work" begin
        renamed = NMRAnalysis.renamedkeyword
        @test renamed([0.1], nothing, :relaxationtimes, :T) == [0.1]
        @test (@test_logs renamed(nothing, [0.1], :relaxationtimes, :T)) == [0.1]
        @test_throws ArgumentError renamed([0.1], [0.2], :relaxationtimes, :T)
        for f in (relaxation2d, recovery2d, methylccr2d, ccr2d, cpmg2d, cest2d, pre2d,
                  lpre2d, relaxation1d, tract, kinetics1d,
                  NMRAnalysis.Analysis1D.nutationexperiment)    # behind calibration1d
            @test any(m -> :T in Base.kwarg_decl(m), methods(f))
        end
    end

    @testset "Which planes are saturated" begin
        @test_throws ArgumentError asksaturation(2; prompt=false)
        # the record of a hetNOE reads back as the saturation it was given
        call = analysiscall("hetnoe2d", ["1", "2"]; saturation=[false, true])
        text = callstring(call)
        @test occursin("saturation=[false, true]", text)
        ex = Meta.parse(text)
        kw = only(ex.args[2].args)
        @test literal(kw.args[2]) == [false, true]
    end
end
