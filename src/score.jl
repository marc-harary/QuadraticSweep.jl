# Score for coefficient of determination (R2)
function score_r2(d::Dataset)
    num = (d.s_xy - 1 / d.n * d.s_x * d.s_y)^2
    den = (d.s_xx - 1 / d.n * d.s_x^2) * (d.s_yy - 1 / d.n * d.s_y^2)
    return num / den
end

# Lift for R2
function lift_r2(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .^ 2, x .* y, y .^ 2, x, y)
end

# Score for Pearson correlation coefficient (cor)
function score_cor(d::Dataset)
    return d.s_xy / sqrt(d.s_xx * d.s_yy + 1e-10)
end

# Life for cor
function lift_cor(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return lift_r2(x, y)
end

# Score for total variation (tv)
function score_tv(d::Dataset)
    return abs(d.s_x - d.s_y) / d.n
end

# Lift for tv
function lift_tv(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .^ 2, y .^ 2, x, y)
end

# Score for difference of variances (dv)
function score_dv(d::Dataset)
    return abs(d.s_xx - d.s_yy) / d.n
end

# Lift for dv 
function lift_dv(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return lift_tv(mat)
end

# Score for covariance (cov)
function score_cov(d::Dataset)
    return (d.s_xy - d.s_x * d.s_y / d.n) / d.n
end

# Lift for cov
function lift_cov(x::Vector{T}, y::Vector{T})::Matrix{T} where {T}
    return hcat(x .* y, x, y)
end

# Mapping from symbols to score functions
const SCORE_FUNCTIONS = Dict(
    :r2 => (score_r2, lift_r2, false),
    :cor => (score_cor, lift_cor, false),
    :tv => (score_tv, lift_tv, true),
    :cov => (score_cov, lift_cov, false),
    :dv => (score_dv, lift_dv, false)
)
