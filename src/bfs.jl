using Random
using Combinatorics
using Statistics

mutable struct Dataset
    s_xx::Float64
    s_xy::Float64
    s_yy::Float64
    s_x::Float64
    s_y::Float64
    n::Int64
    j::Float64
end

function update(ipt::Dataset, p::Tuple{Float64, Float64}, J::Function)
    x, y = p
    opt = Dataset(
        ipt.s_xx + x * x,
        ipt.s_xy + x * y,
        ipt.s_yy + y * y,
        ipt.s_x + x,
        ipt.s_y + y,
        ipt.n + 1,
        0.0
    )
    opt.j = J(opt)
    return opt
end

function bitmask2array(bitmask::Int, num_steps::Int)
    # Initialize an array to store the binary representation
    binary_array = Vector{Int}(undef, num_steps)

    # Extract bits from the bitmask
    for i in 1:num_steps
        # Extract the i-th bit from the right
        binary_array[num_steps - i + 1] = (bitmask >> (i - 1)) & 1
    end

    return Bool.(binary_array)
end

# Define a Node struct for a tree representation of the problem
mutable struct TreeNode
    dataset::Dataset
    inc::Union{TreeNode, Nothing}
    exc::Union{TreeNode, Nothing}
end

# BFS to find the best sequence of includes and excludes using bitmasks
function bfs(data::Array{Float64, 2}, Tc::Array{Int64}, T::Array{Int64}, J::Function)
    # Initialize the root of the tree
    root_dataset = Dataset(0.0, 0.0, 0.0, 0.0, 0.0, 0, -Inf)
    root = TreeNode(root_dataset, nothing, nothing)

    # Queue will store tuples of (current TreeNode, bitmask, Tc, T)
    queue = [(root, 0, Tc, T)]  # bitmask starts at 0

    best_node = root
    best_bitmask = 0

    while !isempty(queue)
        # Pop the first element from the queue
        node, bitmask, current_Tc, current_T = popfirst!(queue)

        # If no more elements in T, stop processing this branch
        if isempty(current_T)
            if node.dataset.j > best_node.dataset.j
                best_node = node
                best_bitmask = bitmask
            end
        else
            # Explore the "include" branch
            if !isempty(current_Tc)
                next_point = Tuple(data[current_Tc[1], :])
                inc_dset = update(node.dataset, next_point, J)
                node.inc = TreeNode(inc_dset, nothing, nothing)
                new_bitmask = (bitmask << 1) | 0  # Set the next bit to 0 (include)
                push!(queue, (node.inc, new_bitmask, current_Tc[2:end], current_T[2:end]))
            end

            # Explore the "exclude" branch
            if !isempty(current_T)
                next_point = Tuple(data[current_T[1], :])
                exc_dset = update(node.dataset, next_point, J)
                node.exc = TreeNode(exc_dset, nothing, nothing)
                new_bitmask = (bitmask << 1) | 1  # Set the next bit to 1 (exclude)
                push!(queue, (node.exc, new_bitmask, current_Tc, current_T[2:end]))
            end
        end
    end

    # Convert bitmask to Boolean array indicating indices used in T
    best_path = bitmask2array(best_bitmask, length(T))
    T_used = T[best_path]

    # Get indices used in T^c
    num_T_used = sum(best_path)
    num_Tc_used = length(T) - num_T_used
    Tc_used = Tc[1:num_Tc_used]

    # Get full set of indices used
    idxs = sort(vcat(T_used, Tc_used))

    return idxs, best_node.dataset.j
end
