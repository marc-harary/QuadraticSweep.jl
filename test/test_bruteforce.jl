using Random
using Statistics
using ProgressMeter
using QuadraticSweep
using Test

# Experiment parameters
seed_rng = MersenneTwister(123)
n_iter = 1_000
k = 10
n = 20
score = :r2

# Wrap the test loop inside a testset
@testset "Quadratic Sweep vs Brute Force" begin
    @showprogress for i in 1:n_iter
        # Seed and generate data
        seed = rand(seed_rng, UInt128)
        data_rng = MersenneTwister(seed)
        x, y = rand(data_rng, n), rand(data_rng, n)

        # Separate programmatically
        idxs_prd, score_prd = sweep(x, y; k = k, score = score)

        # Separate via brute-force
        idxs_grd, score_grd = brute_force(x, y; k = k, score = score)

        # Check for match, log seed if not
        @test idxs_prd == idxs_grd
        @test isapprox(score_prd, score_grd)
    end
end
