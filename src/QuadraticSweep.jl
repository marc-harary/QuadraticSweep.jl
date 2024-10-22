module QuadraticSweep

export sweep, brute_force

export func

"""
    func(x)

Return double the number `x` plus `1`.
"""
func(x) = 2x + 1


include("score.jl")
include("bfs.jl")
include("sweep.jl")

end
