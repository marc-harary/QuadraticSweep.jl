using Combinatorics
using LinearAlgebra
using Ipopt
using JuMP
using Distributed
using Suppressor
using Base.Threads
using PyCall

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

function brute_force(x::Vector{Float64}, y::Vector{Float64}, k::Int,
        score::Function)::Tuple{Vector{Int64}, Vector{Int64}}
    @assert size(x)==size(y) "x and y must be the same shape"
    @assert k<size(x, 1) "k must be less than n"

    n = size(x, 1)
    best_val = -Inf
    best_idxs = nothing

    # Try each subset of size k
    for idxs in combinations(1:n, k)
        if (val = score(x[idxs], y[idxs])) > best_val
            best_val = val
            best_idxs = idxs
        end
    end

    # Get complement of best set
    comp_idxs = setdiff(1:n, best_idxs)

    return best_idxs, comp_idxs
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
