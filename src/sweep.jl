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

    best_score = -Inf
    best_idxs = nothing

    # Try each subset of size k
    for idxs in combinations(1:n, k)
        dataset = Dataset(x[idxs], y[idxs], score_func)
        if dataset.j > best_score
            best_score = dataset.j
            best_idxs = sort(idxs)
        end
    end

    return best_idxs, best_score
end

"""
    sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}

Efficiently identifies the subset of `k` points from the dataset `(x, y)` that maximizes a given score function, using a topological sweep based on a lifting transformation.

# Arguments
- `x::Vector{T}`: A vector of independent variable values, where `T` is a subtype of `Number`.
- `y::Vector{T}`: A vector of dependent variable values, with the same length as `x`.
- `k::Int`: The number of points to select in the optimal subset. Must be between 1 and `n-1` where `n` is the total number of points.
- `score::Symbol`: A symbol representing the score function to use (e.g., `:r2`, `:cor`, `:tv`). The symbol should correspond to a valid function in `SCORE_FUNCTIONS`.

# Returns
- A tuple containing:
    1. `best_idxs::Vector{Int64}`: The indices of the `k` points that form the best subset.
    2. `best_score::Float64`: The maximum score value for the best subset of size `k`.

# Methodology
- The function applies a **lifting transformation** to the data points, projecting them into a higher-dimensional space. This facilitates the identification of an optimal separating hyperplane.
- It then performs a **topological sweep**, iterating over candidate subsets (`k`-sets) to efficiently search for the subset with the highest score.
- If the dataset size `n` is less than or equal to the embedding dimension `d` (from the score function), or if `k <= d`, the function falls back to a brute-force search.

# Notes
- The function uses a combination of **combinatorics** and **linear algebra** to determine the optimal subset.
- The score function and the corresponding lifting transformation are retrieved based on the provided `score` symbol.
- The method ensures efficiency by skipping a brute-force search when possible, instead leveraging the topological sweep to reduce the number of computations.

# Throws
- `AssertionError` if `x` and `y` do not meet the required constraints (`k < n`, `x` and `y` are the same size).
- `ErrorException` if the `score` symbol is invalid or not available in `SCORE_FUNCTIONS`.

# Example
```julia
best_idxs, best_score = sweep([1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 6.0, 8.0]; k = 2, score = :r2)

# See also
- brute_force: The fallback method for an exhaustive search.
- check_input: Validates inputs and retrieves the score function configuration.
"""
function sweep(
        x::Vector{Float64}, y::Vector{Float64}; k::Int64, score::Symbol)::Tuple{
        Vector{Int64}, Float64}
    # Fetch the config
    score_func, lift_func, rev, d, n = check_input(x, y; k, score)

    # Can only do sweep if n > d
    if n <= d || k <= d
        return brute_force(x, y; k = k, score = score)
    end

    # Lift dataset
    lifted = lift_func(x, y)

    best_idxs = nothing
    best_score = -Inf

    # Iterate over k-sets
    for P in combinations(1:n, d)
        # Compute hyperplane
        lifted_P = lifted[P, :]
        lifted_aug = hcat(lifted_P, ones(d))
        N = nullspace(lifted_aug)
        omega = N[1:(end - 1)]

        # Sort points by inner products
        C = setdiff(1:n, P)
        phi = lifted[C, :] * omega
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
