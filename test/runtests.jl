using NMRAnalysis
using Test
using SafeTestsets

@safetestset "diffusion test" begin
    include("diffusion_test.jl")
end

@safetestset "Residue Label Parser Tests" begin
    # # Standard format tests
    # @test NMRAnalysis.GUI2D.extract_residue_number("A7") == 7
    # @test NMRAnalysis.GUI2D.extract_residue_number("D98") == 98
    # @test NMRAnalysis.GUI2D.extract_residue_number("G321") == 321

    # # Reversed format tests
    # @test NMRAnalysis.GUI2D.extract_residue_number("7A") == 7
    # @test NMRAnalysis.GUI2D.extract_residue_number("98D") == 98
    # @test NMRAnalysis.GUI2D.extract_residue_number("321G") == 321

    # # Methyl group tests
    # @test NMRAnalysis.GUI2D.extract_residue_number("I13CD1") == 13
    # @test NMRAnalysis.GUI2D.extract_residue_number("L98CD2") == 98
    # @test NMRAnalysis.GUI2D.extract_residue_number("M98CE") == 98

    # # Non-standard residue tests
    # @test NMRAnalysis.GUI2D.extract_residue_number("X99") == -99
    # @test NMRAnalysis.GUI2D.extract_residue_number("Z100") == -100

    # # Error cases
    # @test_throws ArgumentError NMRAnalysis.GUI2D.extract_residue_number("")
    # @test NMRAnalysis.GUI2D.extract_residue_number("ABC") == 0
    # @test NMRAnalysis.GUI2D.extract_residue_number("N") == 0
    # @test NMRAnalysis.GUI2D.extract_residue_number("X") == 0
end
