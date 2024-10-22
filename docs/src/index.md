# QuadraticSweep.jl Documentation

## Overview

**QuadraticSweep.jl** is a Julia package that provides efficient algorithms for selecting a subset of data points in 2D that maximizes various statistical criteria such as the coefficient of determination ($R^2$), correlation coefficient, total variation, and more. The package leverages combinatorial geometry and optimization techniques to find these subsets.

The package implements two key algorithms:
- **Sweep Algorithm**: Efficient topological sweep that projects data into higher dimensions to identify optimal subsets.
- **Brute Force Algorithm**: Exhaustive search through all possible subsets.

## Installation

To install QuadraticSweep.jl, you can use Julia's package manager. This package is not yet registered in the General registry, so you'll need to install it directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url="https://github.com/username/QuadraticSweep.jl")
```

Ensure that you have Julia installed on your system. If not, download and install it from https://julialang.org/downloads/.

## Main methods
```@docs
sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}
brute_force(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}
```

## Usage 
Once installed, you can start using the QuadraticSweep module. Below is an example of how to use the `sweep` and `brute_force` functions to find the best subset of points maximizing a given score:

```julia
using QuadraticSweep

# Sample data
x = [1.0, 2.0, 3.0, 4.0]
y = [2.0, 3.0, 6.0, 8.0]

# Finding the best subset using the R^2 score
k = 2
score = :r2
best_idxs, best_score = sweep(x, y; k=k, score=score)

println("Best subset indices: ", best_idxs)
println("Best score: ", best_score)
```

## Supported score functions
QuadraticSweep supports several score functions, and you can specify the score function to maximize by passing the appropriate symbol to the `sweep` function:
<table>
    <thead>
        <tr>
            <th>Name</th>
            <th>Symbol</th>
            <th>Score Equation</th>
            <th>Maximizing</th>
            <th>Embedding Equation</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td>Coefficient of Determination</td>
            <td><code>:r2</code></td>
            <td>$R^2 = \frac{\left( S_{XY} - \frac{1}{n} S_X S_Y \right)^2}{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}$</td>
            <td>True</td>
            <td>$\mathcal L(x, y) = \left( x^2, xy, y^2, x, y \right)$</td>
        </tr>
        <tr>
            <td>Correlation Coefficient</td>
            <td><code>:cor</code></td>
            <td>$r = \frac{S_{XY} - \frac{1}{n} S_X S_Y}{\sqrt{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}}$</td>
            <td>True</td>
            <td>$\mathcal L(x, y) = \left( x^2, xy, y^2, x, y \right)$</td>
        </tr>
        <tr>
            <td>Total Variation</td>
            <td><code>:tv</code></td>
            <td>$TV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) + \left( S_{YY} - \frac{1}{k} S_Y^2 \right)$</td>
            <td>False</td>
            <td>$\mathcal L(x, y) = \left( x^2, y^2, x, y \right)$</td>
        </tr>
        <tr>
            <td>Covariance</td>
            <td><code>:cov</code></td>
            <td>$\text{cov}(X, Y) = S_{XY} - \frac{1}{k} S_X S_Y$</td>
            <td>True</td>
            <td>$\mathcal L(x, y) = \left( x, y, xy \right)$</td>
        </tr>
        <tr>
            <td>Difference of Variances</td>
            <td><code>:dv</code></td>
            <td>$DV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) - \left( S_{YY} - \frac{1}{k} S_Y^2 \right)$</td>
            <td>True</td>
            <td>$\mathcal L(x, y) = \left( x^2, y^2, x, y \right)$</td>
        </tr>
        <tr>
            <td>Fraction of Variance Unexplained</td>
            <td><code>:fvu</code></td>
            <td>$FVU = \frac{S_{YY}}{\left( S_{XX} - \frac{1}{k} S_X^2 \right)}$</td>
            <td>False</td>
            <td>$\mathcal L(x, y) = \left( x, x^2, y^2 \right)$</td>
        </tr>
    </tbody>
</table>
