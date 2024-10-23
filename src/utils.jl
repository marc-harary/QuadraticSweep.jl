function check_input(x::Vector{T}, y::Vector{T}; k::Int,
        score::Symbol)::Tuple{Function, Function, Bool, Int64, Int64} where {T <: Number}
    n = length(x)
    @assert length(y)==n "x and y must be the same shape"
    @assert k<n "k must be less than n"
    config = get(SCORE_FUNCTIONS, score, nothing)
    if config === nothing
        error("Invalid score symbol: $score. Available options are: $(keys(SCORE_FUNCTIONS))")
    end
    return config..., n
end

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
