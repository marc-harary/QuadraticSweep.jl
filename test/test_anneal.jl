using Random
using Statistics
using ProgressMeter
using Suppressor

include("qp.jl")

# Experiment parameters
init_seed = 1234
seed_rng = MersenneTwister(init_seed)
steps = 10_000
Tmax = 100_000
Tmin = 1
n_iter = 1_000
k = 13
n = 25
score = (x, y) -> cor(x, y)^2

matches = []
ratios = []
@showprogress for i in 1:n_iter
    # Seed and generate data
    seed = rand(seed_rng, UInt128)
    data_rng = MersenneTwister(seed)
    x, y = rand(data_rng, n), rand(data_rng, n)

    # Separate programatically
    idxs_prd = find_optimal_subset(
        x, y, k, score; steps = steps, seed = seed, Tmax = Tmax, Tmin = Tmin)
    score_prd = score(x[idxs_prd], y[idxs_prd])

    # Separate via brute-force
    idxs_grd, _ = brute_force(x, y, k, score)
    score_grd = score(x[idxs_grd], y[idxs_grd])

    push!(ratios, score_prd / score_grd)
    push!(matches, sort(idxs_prd) == sort(idxs_grd))
end

println("Average success rate: ", mean(matches))
println("Average R2 ratio: ", mean(ratios))
