using Random
using Statistics
using ProgressMeter
using QuadraticSweep
using Test

# Experiment parameters
init_seed = 123
seed_rng = MersenneTwister(init_seed)
n_iter = 1_000
k = 10
n = 20
lift = (x, y) -> hcat(x, y, x .^ 2, y .^ 2, x .* y)
score = (x, y) -> cor(x, y)^2
function score_struct(d)
    num = (d.s_xy - 1 / d.n * d.s_x * d.s_y)^2
    den = (d.s_xx - 1 / d.n * d.s_x^2) * (d.s_yy - 1 / d.n * d.s_y^2)
    return num / den
end

# Wrap the test loop inside a testset
@testset "Quadratic Sweep vs Brute Force" begin
    @showprogress for i in 1:n_iter
        # Seed and generate data
        seed = rand(seed_rng, UInt128)
        data_rng = MersenneTwister(seed)
        x, y = rand(data_rng, n), rand(data_rng, n)

        # Separate programmatically
        # idxs_prd = sweep(x, y; k = k, L = lift, S = score_struct, rho = false)
        idxs_prd = sweep(x, y; k = k, score = :r2)

        # Separate via brute-force
        idxs_grd, _ = brute_force(x, y, k, score)
        r2 = cor(x[idxs_grd], y[idxs_grd])^2

        # Check for match, log seed if not
        @test sort(idxs_prd) == sort(idxs_grd)
    end
end
