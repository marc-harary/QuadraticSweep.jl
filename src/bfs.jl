# Convert a bitmask to a Boolean array
function bitmask2array(bitmask::Int, num_steps::Int)
    # Initialize an array to store the binary representation
    binary_array = Vector{Bool}(undef, num_steps)

    # Extract bits from the bitmask
    for i in 1:num_steps
        # Extract the i-th bit from the right
        binary_array[num_steps - i + 1] = Bool((bitmask >> (i - 1)) & 1)
    end

    return binary_array
end

# TreeNode structure
mutable struct TreeNode
    dataset::Dataset
    inc::Union{TreeNode, Nothing}
    exc::Union{TreeNode, Nothing}
end

# Performs breadth-first search (BFS) to get optimal set of inliers
function bfs(
        data::Array{Float64, 2}, k::Int64, Tc::Array{Int64}, T::Array{Int64},
        J::Function)::Tuple{Vector{Int64}, Float64}
    d = length(T)  # Length of the exclude set T

    # Initialize root dataset with the first k-d points in data
    x = data[Tc[1:(k - d)], 1]
    y = data[Tc[1:(k - d)], 2]
    root_dataset = Dataset(x, y, J)

    # Update current T^c so we don't double-count first k-d points
    Tc_og = copy(Tc)
    Tc = Tc[(k - d + 1):end]

    # Create the root node
    root = TreeNode(root_dataset, nothing, nothing)

    # Initialize the queue with the root node, starting bitmask, and current Tc and T
    queue = [(root, 0, Tc, T)]

    # Best argument and score
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
                new_bitmask = (bitmask << 1) | 0
                push!(queue, (node.inc, new_bitmask, current_Tc[2:end], current_T[2:end]))
            end

            # Explore the "exclude" branch
            if !isempty(current_T)
                next_point = Tuple(data[current_T[1], :])
                exc_dset = update(node.dataset, next_point, J)
                node.exc = TreeNode(exc_dset, nothing, nothing)
                new_bitmask = (bitmask << 1) | 1
                push!(queue, (node.exc, new_bitmask, current_Tc, current_T[2:end]))
            end
        end
    end

    # Reconstruct original set of indices used
    best_path = bitmask2array(best_bitmask, d)
    T_used = T[best_path]
    Tc_used = Tc_og[1:(k - sum(best_path))]
    best_idxs = sort(vcat(T_used, Tc_used))
    best_score = best_node.dataset.j

    return best_idxs, best_score
end
