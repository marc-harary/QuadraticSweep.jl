using Test
using Random
using QuadraticSweep

include("test_cases.jl")

# Unit tests for error handling and basic functionality
@testset "Quadratic sweep tests" begin

    # Test input length mismatch error
    @testset "Error-handling" begin
        x = [1.0, 2.0, 3.0]
        y = [4.0, 5.0]
        k = 2

        # x and y length mismatch
        @test_throws AssertionError sweep(x, y; k = k, score = :r2)

        # Invalid k (k > n)
        y = [4.0, 5.0, 6.0]
        @test_throws AssertionError sweep(x, y; k = 4, score = :dv)

        # Invalid score symbol
        @test_throws ErrorException sweep(x, y; k = 2, score = :invalid_symbol)
    end

    @testset "Sweep correctness" begin
        for (score, n, k, seed, idxs_grd, val_grd) in test_cases
            # Generate data from random seed
            data_rng = MersenneTwister(seed)
            x, y = rand(data_rng, n), rand(data_rng, n)

            # Run algorithm and check indices and score
            idxs_sweep, val_sweep = sweep(x, y; k = k, score = score)
            @test idxs_grd == idxs_sweep
            @test isapprox(val_grd, val_sweep)
        end
    end

    @testset "Brute-force correctness" begin
        for (score, n, k, seed, idxs_grd, val_grd) in test_cases
            # Generate data from random seed
            data_rng = MersenneTwister(seed)
            x, y = rand(data_rng, n), rand(data_rng, n)

            # Run algorithm and check indices and score 
            idxs_brute, val_brute = brute_force(x, y; k = k, score = score)
            @test idxs_grd == idxs_brute
            @test isapprox(val_grd, val_brute)
        end
    end
end
