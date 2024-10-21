using Combinatorics
using LinearAlgebra

function sweep(
        X::Vector, Y::Vector; k::Int64, L::Function, S::Function, rho::Bool)::Vector{Int64}
    # Lift dataset and store embedding dimension
    LX = L(X, Y)
    n, d = size(LX)

    alpha_max = -Inf
    S_max = nothing

    # Iterate over k-sets
    for P in combinations(1:n, d)
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

        # Do BFS on current set
        S_cur, alpha_cur = bfs(hcat(X, Y), k, C, P, S)

        if alpha_cur > alpha_max
            alpha_max = alpha_cur
            S_max = S_cur
        end
    end

    return S_max
end
