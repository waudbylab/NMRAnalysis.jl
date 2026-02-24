using NMRAnalysis
using NMRAnalysis: SingleElementVector, StandardVector
using Test

@testset "MaybeVector" begin
    @testset "Construction" begin
        s = MaybeVector(1)
        @test s isa SingleElementVector{Int}
        @test s.x == 1

        v = MaybeVector([1, 2, 3])
        @test v isa StandardVector{Int}
        @test v.x == [1, 2, 3]

        sf = MaybeVector(1.0)
        @test sf isa SingleElementVector{Float64}

        vf = MaybeVector([1.0, 2.0])
        @test vf isa StandardVector{Float64}
    end

    @testset "Indexing" begin
        s = MaybeVector(1)
        @test s[1] == 1
        @test s[2] == 1  # All valid indices return the same value

        v = MaybeVector([1, 2, 3])
        @test v[1] == 1
        @test v[2] == 2
        @test v[3] == 3
        @test_throws BoundsError v[4]
    end

    @testset "Length and Size" begin
        s = MaybeVector(1)
        @test length(s) == 1
        @test size(s) == (1,)

        v = MaybeVector([1, 2, 3])
        @test length(v) == 3
        @test size(v) == (3,)
    end

    @testset "Iteration" begin
        s = MaybeVector(1)
        collected_s = collect(s)
        @test collected_s == [1]

        v = MaybeVector([1, 2, 3])
        collected_v = collect(v)
        @test collected_v == [1, 2, 3]
    end

    @testset "Broadcasting" begin
        s = MaybeVector(2)
        @test s .* 2 == 4

        v = MaybeVector([1, 2, 3])
        @test v .* 2 == [2, 4, 6]
    end

    @testset "Conversion" begin
        # Single -> Standard
        s = MaybeVector(1)
        converted_to_standard = convert(StandardVector{Int}, s)
        @test converted_to_standard isa StandardVector{Int}
        @test converted_to_standard.x == [1]

        # Standard -> Single (valid case)
        v = MaybeVector([1])
        converted_to_single = convert(SingleElementVector{Int}, v)
        @test converted_to_single isa SingleElementVector{Int}
        @test converted_to_single.x == 1

        # Standard -> Single (invalid case)
        v_long = MaybeVector([1, 2])
        @test_throws ArgumentError convert(SingleElementVector{Int}, v_long)

        # To Vector
        @test convert(Vector{Int}, s) == [1]
        @test convert(Vector{Int}, v) == [1]

        # From Vector
        @test convert(MaybeVector{Int}, [1]) isa SingleElementVector{Int}
        @test convert(MaybeVector{Int}, [1, 2]) isa StandardVector{Int}
    end

    @testset "Mutation" begin
        s = MaybeVector(1)
        s[1] = 2
        @test s[1] == 2
        @test s[2] == 2  # All indices show the new value

        v = MaybeVector([1, 2, 3])
        v[2] = 4
        @test v[2] == 4
        @test v.x == [1, 4, 3]
    end

    @testset "Copy and Similar" begin
        s = MaybeVector(42)
        sc = copy(s)
        @test sc isa SingleElementVector{Int}
        @test sc.x == 42
        sc[1] = 99
        @test s[1] == 42  # original unchanged

        v = MaybeVector([1, 2, 3])
        vc = copy(v)
        @test vc isa StandardVector{Int}
        @test vc.x == [1, 2, 3]
        vc[1] = 99
        @test v[1] == 1  # original unchanged

        ss = similar(s)
        @test ss isa SingleElementVector{Int}

        vs = similar(v)
        @test vs isa StandardVector{Int}
        @test length(vs) == 3
    end

    @testset "IndexStyle" begin
        @test Base.IndexStyle(SingleElementVector{Int}) == IndexLinear()
        @test Base.IndexStyle(StandardVector{Int}) == IndexLinear()
    end
end
