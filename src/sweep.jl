using Combinatorics
using LinearAlgebra

function check_input(x::Vector{Float64}, y::Vector{Float64}; k::Int64,
        score::Symbol)::Tuple{Function, Function, Bool, Int64, Int64}
    n = length(x)
    @assert length(y)==n "x and y must be the same shape"
    @assert k<n "k must be less than n"
    config = get(SCORE_FUNCTIONS, score, nothing)
    if config === nothing
        error("Invalid score symbol: $score. Available options are: $(keys(SCORE_FUNCTIONS))")
    end
    return config..., n
end

function brute_force(x::Vector{Float64}, y::Vector{Float64}; k::Int64,
        score::Symbol)::Tuple{Vector{Int64}, Float64}
    # Fetch the config
    score_func, _, _, _, n = check_input(x, y; k, score)

    best_val = -Inf
    best_idxs = nothing

    # Try each subset of size k
    for idxs in combinations(1:n, k)
        dataset = Dataset(x[idxs], y[idxs], score_func)
        if dataset.j > best_val
            best_val = dataset.j
            best_idxs = sort(idxs)
        end
    end

    return best_idxs, best_val
end

function sweep(
        x::Vector{Float64}, y::Vector{Float64}; k::Int64, score::Symbol)::Tuple{
        Vector{Int64}, Float64}
    # Fetch the config
    score_func, lift_func, rev, d, n = check_input(x, y; k, score)

    # Can only do sweep if n > d
    if n <= d
        return brute_force(x, y; k = k, score = score)
    end

    # Lift dataset
    lifted = lift_func(x, y)

    best_idxs = -Inf
    best_score = nothing

    # Iterate over k-sets
    for P in combinations(1:n, d)
        # Compute hyperplane
        lifted_P = lifted[P, :]
        lifted_aug = hcat(lifted_P, ones(d))
        N = nullspace(lifted_aug)
        omega = N[1:(end - 1)]

        # Sort points by inner products
        C = setdiff(1:n, P)
        phi = LX[C, :] * omega
        I = sortperm(phi, rev = rev)
        C = C[I]

        # Do BFS on current set
        cur_idxs, cur_score = bfs(hcat(x, y), k, C, P, score_func)

        # Update best values if necessary
        if cur_score > best_score
            best_score = cur_score
            best_idxs = cur_idxs
        end
    end

    return best_idxs, best_score
end
