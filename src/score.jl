# Constant to prevent division by zero
const EPS = 1e-16

# Dataset structure
mutable struct Dataset
    s_xx::Float64
    s_xy::Float64
    s_yy::Float64
    s_x::Float64
    s_y::Float64
    n::Int64
    j::Float64
    function Dataset(
            x::Vector{Float64}, y::Vector{Float64}, J::Union{Function, Nothing} = nothing)
        opt = new(sum(x .^ 2), sum(x .* y), sum(y .^ 2), sum(x), sum(y), length(x), 0)
        if !isnothing(J)
            opt.j = J(opt)
        end
        return opt
    end
    function Dataset(s_xx::Float64, s_xy::Float64, s_yy::Float64, s_x::Float64,
            s_y::Float64, n::Int64, J::Union{Function, Nothing} = nothing)
        opt = new(s_xx, s_xy, s_yy, s_x, s_y, n, 0)
        if !isnothing(J)
            opt.j = J(opt)
        end
        return opt
    end
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
        J
    )
    return opt
end

# Score for coefficient of determination (R2)
function score_r2(d::Dataset)
    num = (d.s_xy - 1 / d.n * d.s_x * d.s_y)^2
    den = (d.s_xx - 1 / d.n * d.s_x^2) * (d.s_yy - 1 / d.n * d.s_y^2)
    return num / (den + EPS)
end

# Lift for R2
function lift_r2(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .^ 2, x .* y, y .^ 2, x, y)
end

# Score for Pearson correlation coefficient (cor)
function score_cor(d::Dataset)
    num = d.s_xy - 1 / d.n * d.s_x * d.s_y
    den = sqrt(d.s_xx - 1 / d.n * d.s_x^2) * sqrt(d.s_yy - 1 / d.n * d.s_y^2)
    return num / (den + EPS)
end

# Life for cor
function lift_cor(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return lift_r2(x, y)
end

# Score for total variation (tv)
function score_tv(d::Dataset)
    var_x = d.s_xx / d.n - (d.s_x / d.n)^2
    var_y = d.s_yy / d.n - (d.s_y / d.n)^2
    return var_x + var_y
end

# Lift for tv
function lift_tv(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .^ 2, y .^ 2, x, y)
end

# Score for difference of variances (dv)
function score_dv(d::Dataset)
    var_x = d.s_xx / d.n - (d.s_x / d.n)^2
    var_y = d.s_yy / d.n - (d.s_y / d.n)^2
    return var_x - var_y
end

# Lift for dv 
function lift_dv(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return lift_tv(x, y)
end

# Score for covariance (cov)
function score_cov(d::Dataset)
    return d.s_xy / d.n - d.s_x * d.s_y / d.n^2
end

# Lift for cov
function lift_cov(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .* y, x, y)
end

# Score for fraction of variance unexplained (fvu)
function score_fvu(d::Dataset)
    return d.s_yy / (d.s_xx - 1 / d.n * d.s_x^2)
end

# Lift for fvu
function lift_fvu(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x, x .^ 2, y .^ 2)
end

# Mapping from symbols to score functions
const SCORE_FUNCTIONS = Dict(
    :r2 => (score_r2, lift_r2, false, 5),
    :cor => (score_cor, lift_cor, false, 5),
    :tv => (score_tv, lift_tv, true, 4),
    :cov => (score_cov, lift_cov, false, 3),
    :dv => (score_dv, lift_dv, false, 4),
    :fvu => (score_fvu, lift_fvu, true, 3)
)
