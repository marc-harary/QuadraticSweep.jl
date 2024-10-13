using Random
using Statistics
using ProgressMeter
using DataFrames
using QuadraticSweep
import MathOptInterface as MOI

# Experiment parameters
tol = 1e-10
init_seed = 1234
seed_rng = MersenneTwister(init_seed)
n_iter = 1_000
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

scores = [
    (x, y) -> cor(x, y),
    (x, y) -> cor(x, y)^2,
    (x, y) -> cov(x),
    (x, y) -> var(x) + var(y),
    (x, y) -> var(x) - var(y),
]

lifts = [
    (x, y) -> hcat(x.^2, x.*y, y.^2, x, y),
    (x, y) -> hcat(x.^2, y.^2, x, y),
    (x, y) -> hcat(x.*y, x, y),
]

# DataFrame to store the means and standard deviations of all combinations
summary_df = DataFrame(
    score_func = String[],
    lift_func = String[],
    mean_success = Float64[],
    std_success = Float64[],
    mean_objective_value = Float64[],
    std_objective_value = Float64[],
    mean_dual_objective_value = Float64[],
    std_dual_objective_value = Float64[],
    mean_barrier_iterations = Float64[],
    std_barrier_iterations = Float64[],
    mean_solve_time = Float64[],
    std_solve_time = Float64[]
)

for score_idx in 1:length(scores)
    for lift_idx in 1:length(lifts)
        score = scores[score_idx]
        lift = lifts[lift_idx]

        results_df = DataFrame(
            success_value = Bool[],
            objective_value = Float64[],
            dual_objective_value = Float64[],
            barrier_iterations = Int[],
            solve_time = Float64[],
        )

        @showprogress for i=1:n_iter
            # Seed and generate data
            seed = rand(seed_rng, UInt128)
            data_rng = MersenneTwister(seed)
            x, y = rand(data_rng, n), rand(data_rng, n)

            # Separate via brute-force
            best_idxs, comp_idxs = brute_force(x, y, k, score)

            # Lift data
            lifted = lift(x, y)
            A = lifted[best_idxs, :]
            B = lifted[comp_idxs, :]

            model = lp_intersect(A, B)

            # Collect optimization attributes
            obj_value = MOI.get(model, MOI.ObjectiveValue())
            dual_obj_value = MOI.get(model, MOI.DualObjectiveValue())
            barrier_iters = MOI.get(model, MOI.BarrierIterations())
            solve_time = MOI.get(model, MOI.SolveTimeSec())
            success = obj_value > tol

            # Append results to the DataFrame
            push!(results_df, (success, obj_value, dual_obj_value, barrier_iters, solve_time))
        end

        # Compute mean and standard deviation of the attributes
        mean_values = combine(results_df, names(results_df) .=> mean)
        std_values = combine(results_df, names(results_df) .=> std)

        # Append the mean and std results to the summary DataFrame
        push!(summary_df, (
            "score_$score_idx",  # Adding labels for each score and lift
            "lift_$lift_idx",
            mean_values[1, :success_value_mean], std_values[1, :success_value_std],
            mean_values[1, :objective_value_mean], std_values[1, :objective_value_std],
            mean_values[1, :dual_objective_value_mean], std_values[1, :dual_objective_value_std],
            mean_values[1, :barrier_iterations_mean], std_values[1, :barrier_iterations_std],
            mean_values[1, :solve_time_mean], std_values[1, :solve_time_std]
        ))
    end
end

# Print the summary DataFrame containing the mean and standard deviation of all combinations
println(summary_df)
