using Combinatorics
using LinearAlgebra
using Ipopt
using JuMP
using Distributed
using Suppressor
using Base.Threads
# using Metaheuristics
using PyCall

pushfirst!(pyimport("sys")."path", "");
subset_anneal = pyimport("subset_anneal")

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


############### ORIGINAL ############### 
# function lin_sep(x::Vector, y::Vector; k::Int64, lift::Function, score::Function, rev::Bool)::Vector{Int64}
#     # Lift dataset and store embeddin dimension
#     lifted = lift(x, y)
#     n, m = size(lifted)
# 
#     max_score = -Inf
#     max_set = nothing
# 
#     # Iterate over k-sets
#     for comb in combinations(1:n, m)
#         # Compute hyperplane
#         subset = lifted[comb, :]
#         subset_aug = hcat(subset, ones(m))
#         null_space = nullspace(subset_aug)
#         weights = null_space[1:(end - 1)]
# 
#         # Sort points by inner products
#         comp = setdiff(1:n, comb)
#         prods = lifted[comp, :] * weights
#         idxs = sortperm(prods, rev = rev)
#         idxs = comp[idxs]
# 
#         # Build candidate subsets
#         idx_subsets = [idxs[1:k]]
#         for i in 1:m
#             for comb2 in combinations(comb, i)
#                 push!(idx_subsets, vcat(idxs[1:(k - i)], comb2))
#             end
#         end
# 
#         # Iterate over all subsets
#         for idx_subset in idx_subsets
#             cur_score = score(x[idx_subset], y[idx_subset])
#             if cur_score > max_score
#                 max_score = cur_score
#                 max_set = idx_subset
#             end
#         end
#     end
# 
#     return max_set
# end
######################################## 


function lin_sep(X::Vector, Y::Vector; k::Int64, L::Function, S::Function, rho::Bool)::Vector{Int64}
    # Lift dataset and store embedding dimension
    LX = L(X, Y)
    n, d = size(LX)
    
    alpha_max = -Inf
    S_max = nothing
    
    # Iterate over k-sets
    for P in collect(combinations(1:n, d))
        # Compute hyperplane
        LX_P = LX[P, :]
        LX_aug = hcat(LX_P, ones(d))
        N = nullspace(LX_aug)
        omega = N[1:(end - 1)]
        
        # Sort points by inner products
        C = setdiff(1:n, P)
        phi = LX[C, :] * omega
        I = sortperm(phi, rev = rho)
        C = C[I]
        
        # Iterate over candidate subsets
        for i in 1:Int(ceil(d / 2))
            for P_prime in combinations(P, i)
                S_prime = vcat(C[1:(k - i)], P_prime)
                alpha_cur = S(X[S_prime], Y[S_prime])
                if alpha_cur > alpha_max
                    alpha_max = alpha_cur
                    S_max = S_prime
                end
            end
        end
    end
    
    return S_max
end

function intersect(A::Matrix{Float64}, B::Matrix{Float64}, tol::Float64 = 1e-10)
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



