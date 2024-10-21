using Test
using QuadraticSweep

# Unit tests for error handling and basic functionality
@testset "Sweep and Brute Force Tests" begin

    # Test input length mismatch error
    @testset "check_input Tests" begin
        x = [1.0, 2.0, 3.0]
        y = [4.0, 5.0]
        k = 2

        # x and y length mismatch
        @test_throws AssertionError sweep(x, y; k = k, score = :r2)

        # Invalid k (k >= n)
        y = [4.0, 5.0, 6.0]
        @test_throws AssertionError sweep(x, y; k = 4, score = :r2)

        # Invalid score symbol
        @test_throws ErrorException sweep(x, y; k = 2, score = :invalid_symbol)
    end

    # Test brute force algorithm
    @testset "brute_force Tests" begin
        x = [1.0, 2.0, 3.0]
        y = [4.0, 5.0, 6.0]
        k = 2

        # Simple brute force test
        idxs, score = brute_force(x, y; k, score = :tv)
        @test length(idxs) == k
    end

    # Test sweep algorithm
    @testset "sweep Tests" begin
        x = [1.0, 2.0, 3.0, 4.0]
        y = [2.0, 3.0, 4.0, 5.0]
        k = 3

        # Test sweep for valid inputs
        idxs, score = sweep(x, y; k, score = :dv)
        @test length(idxs) == k
    end
end
