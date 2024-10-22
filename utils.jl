using Combinatorics
using Ipopt
using JuMP
using Distributed
using PyCall
using Suppressor

function find_optimal_subset(x, y, k, score_function; seed, steps, Tmax, Tmin)
    @suppress begin
        optimizer = subset_anneal.SubsetOptimization(x, y, k, score_function, seed)
        optimizer.steps = steps
        optimizer.Tmax = Tmax
        optimizer.Tmin = Tmin
        best_idxs, _ = optimizer.anneal()
        best_idxs .+= 1
        return best_idxs
    end
end

function lp_intersect(A::Matrix{Float64}, B::Matrix{Float64}, tol::Float64 = 1e-10)
    @suppress begin
        @assert size(A, 2)==size(B, 2) "Points must have same dimension"

        m, d = size(A)
        n, _ = size(B)

        # Create model
        model = Model(Ipopt.Optimizer)
        set_silent(model)

        # Define variables
        @variable(model, lambda[1:m]>=0)
        @variable(model, mu[1:n]>=0)

        # Constraints
        @constraint(model, sum(lambda)==1)
        @constraint(model, sum(mu)==1)

        # Objective function
        diff = lambda' * A .- mu' * B
        @objective(model, Min, sum(diff .^ 2))

        # Optimize the model
        optimize!(model)

        # Return whether less than tol
        return model
    end
end

