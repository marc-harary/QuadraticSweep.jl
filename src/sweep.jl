using Combinatorics
using LinearAlgebra
using Ipopt
using JuMP
using Distributed
using Suppressor
using Base.Threads
using PyCall

function sweep(X::Vector, Y::Vector; k::Int64, L::Function, S::Function, rho::Bool)::Vector{Int64}
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
