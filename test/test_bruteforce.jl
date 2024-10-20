using Random
using Statistics
using ProgressMeter
using QuadraticSweep
using Test

# Experiment parameters
init_seed = 1234
seed_rng = MersenneTwister(init_seed)
n_iter = 10_000
k = 8
n = 15
score = (x, y) -> cor(x, y)^2
lift = (x, y) -> hcat(x, y, x .^ 2, y .^ 2, x .* y)

# Convenience function for logging seeds at which experiment failed
function log_failure(seed::UInt128)
    open("seed_log.txt", "a") do file
        write(file, string(seed) * "\n")
    end
end

# Wrap the test loop inside a testset
@testset "Quadratic Sweep vs Brute Force" begin
    @showprogress for i in 1:n_iter
        # Seed and generate data
        seed = rand(seed_rng, UInt128)
        data_rng = MersenneTwister(seed)
        x, y = rand(data_rng, n), rand(data_rng, n)

        # Separate programmatically
        idxs_prd = sweep(x, y; k = k, L = lift, S = score, rho = false)

        # Separate via brute-force
        idxs_grd, _ = brute_force(x, y, k, score)

        # Check for match, log seed if not
        try
            @test sort(idxs_prd) == sort(idxs_grd)
        catch e
            # Log the failure if the test fails
            log_failure(seed)
            @error("Test failed at seed: $seed")
        end
    end
end