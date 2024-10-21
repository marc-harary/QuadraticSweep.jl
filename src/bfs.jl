using Random
using Combinatorics
using Statistics

# Dataset structure
mutable struct Dataset
    s_xx::Float64
    s_xy::Float64
    s_yy::Float64
    s_x::Float64
    s_y::Float64
    n::Int64
    j::Float64
end

# Function to update the dataset with a new point (x, y)
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

# Convert a bitmask to a Boolean array
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

# TreeNode structure
mutable struct TreeNode
    dataset::Dataset
    inc::Union{TreeNode, Nothing}
    exc::Union{TreeNode, Nothing}
end

function bfs(
        data::Array{Float64, 2}, k::Int64, Tc::Array{Int64}, T::Array{Int64}, J::Function)
    d = length(T)  # Length of the exclude set T

    # Initialize root dataset with the first k-d points in data
    x = data[Tc[1:(k - d)], 1]
    y = data[Tc[1:(k - d)], 2]
    root_dataset = Dataset(sum(x .^ 2), sum(x .* y), sum(y .^ 2), sum(x), sum(y), k - d, 0)
    root_dataset.j = J(root_dataset)

    # Update current T^c so we don't double-count first k-d points
    Tc_og = copy(Tc)
    Tc = Tc[(k - d + 1):end]

    # Create the root node
    root = TreeNode(root_dataset, nothing, nothing)

    # Initialize the queue with the root node, starting bitmask, and current Tc and T
    queue = [(root, [], Tc, T)]

    best_node = root
    best_bitmask = nothing
    best_score = -Inf

    while !isempty(queue)
        node, bitmask, current_Tc, current_T = popfirst!(queue)

        # If no more elements in T, update the best node
        if isempty(current_T)
            if node.dataset.j > best_score
                best_node = node
                best_bitmask = copy(bitmask)
                best_score = node.dataset.j
            end
        else
            # Explore the "include" branch
            if !isempty(current_Tc)
                next_point = Tuple(data[current_Tc[1], :])
                inc_dset = update(node.dataset, next_point, J)
                node.inc = TreeNode(inc_dset, nothing, nothing)
                new_bitmask = copy(bitmask)
                push!(new_bitmask, current_Tc[1])
                push!(queue, (node.inc, new_bitmask, current_Tc[2:end], current_T[2:end]))
            end

            # Explore the "exclude" branch
            if !isempty(current_T)
                next_point = Tuple(data[current_T[1], :])
                exc_dset = update(node.dataset, next_point, J)
                node.exc = TreeNode(exc_dset, nothing, nothing)
                new_bitmask = copy(bitmask)
                push!(new_bitmask, current_T[1])
                push!(queue, (node.exc, new_bitmask, current_Tc, current_T[2:end]))
            end
        end
    end

    idxs = vcat(Tc_og[1:(k - d)], best_bitmask)

    return sort(idxs), best_node.dataset.j
end
