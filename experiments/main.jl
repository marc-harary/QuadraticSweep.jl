using Random
using Statistics
using ProgressMeter

include("qp.jl")

# Experiment parameters
init_seed = 1234
seed_rng = MersenneTwister(init_seed)
n_iter = 100_000
k = 8
n = 15
score = (x, y) -> cor(x, y)
lift = (x, y) -> hcat(x, y, x.^2, y.^2, x .* y)

# Convenience function for logging seeds at which experiment failed
function log_failure(seed::UInt128)
    open("seed_log.txt", "a") do file
        write(file, string(seed) * "\n")
    end
end

@showprogress for i=1:n_iter
    # Seed and generate data
    seed = rand(seed_rng, UInt128)
    data_rng = MersenneTwister(seed)
    x, y = rand(data_rng, n), rand(data_rng, n)

    # Separate programatically
    idxs_prd = lin_sep(x, y; k=k, L=lift, S=score, rho=false)
    
    # Separate via brute-force
    idxs_grd, _ = brute_force(x, y, k, score)
    
    # Check for match and log seed if not
    if sort(idxs_prd) != sort(idxs_grd)
        println(seed)
        log_failure(seed)
    end
end
