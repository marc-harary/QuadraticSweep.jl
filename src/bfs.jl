# TreeNode structure
mutable struct TreeNode
    dataset::Dataset
    inc::Union{TreeNode,Nothing}
    exc::Union{TreeNode,Nothing}
end

# Performs breadth-first search (BFS) to get optimal set of inliers
function bfs(
    data::Array{Float64,2},
    k::Int64,
    comp_idxs::Array{Int64},
    pivot_idxs::Array{Int64},
    score_func::Function,
)::Tuple{Vector{Int64},Float64}
    d = length(pivot_idxs)  # Length of the exclude set pivot_idxs

    # Initialize root dataset with the first k-d points in data
    x = data[comp_idxs[1:(k-d)], 1]
    y = data[comp_idxs[1:(k-d)], 2]
    root_dataset = Dataset(x, y, score_func)

    # Update current pivot_idxs so we don't double-count first k-d points
    comp_idxs_og = copy(comp_idxs)
    comp_idxs = comp_idxs[(k-d+1):end]

    # Create the root node
    root = TreeNode(root_dataset, nothing, nothing)

    # Initialize the queue with the root node, starting bitmask, and current comp_idxs and pivot_idxs
    queue = [(root, 0, comp_idxs, pivot_idxs)]

    # Best argument and score
    best_node = root
    best_bitmask = nothing
    best_score = -Inf

    while !isempty(queue)
        node, bitmask, current_comp_idxs, current_pivot_idxs = popfirst!(queue)

        # If no more elements in pivot_idxs, update the best node
        if isempty(current_pivot_idxs)
            if node.dataset.j > best_score
                best_node = node
                best_bitmask = copy(bitmask)
                best_score = node.dataset.j
            end
        else
            # Explore the "include" branch
            if !isempty(current_comp_idxs)
                next_point = Tuple(data[current_comp_idxs[1], :])
                inc_dset = update(node.dataset, next_point, score_func)
                node.inc = TreeNode(inc_dset, nothing, nothing)
                new_bitmask = (bitmask << 1) | 0
                push!(
                    queue,
                    (
                        node.inc,
                        new_bitmask,
                        current_comp_idxs[2:end],
                        current_pivot_idxs[2:end],
                    ),
                )
            end

            # Explore the "exclude" branch
            if !isempty(current_pivot_idxs)
                next_point = Tuple(data[current_pivot_idxs[1], :])
                exc_dset = update(node.dataset, next_point, score_func)
                node.exc = TreeNode(exc_dset, nothing, nothing)
                new_bitmask = (bitmask << 1) | 1
                push!(
                    queue,
                    (node.exc, new_bitmask, current_comp_idxs, current_pivot_idxs[2:end]),
                )
            end
        end
    end

    # Reconstruct original set of indices used
    best_path = bitmask2array(best_bitmask, d)
    pivot_idxs_used = pivot_idxs[best_path]
    comp_idxs_used = comp_idxs_og[1:(k-sum(best_path))]
    best_idxs = sort(vcat(pivot_idxs_used, comp_idxs_used))
    best_score = best_node.dataset.j

    return best_idxs, best_score
end
