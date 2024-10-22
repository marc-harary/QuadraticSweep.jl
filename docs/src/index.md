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
Pkg.add(url="https://github.com/marc-harary/QuadraticSweep.jl")
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

| Name                            | Symbol | Score Equation                                                                                                                                      | Maximizing | Embedding Equation                                       |
|----------------------------------|--------|-----------------------------------------------------------------------------------------------------------------------------------------------------|------------|----------------------------------------------------------|
| Coefficient of Determination     | `:r2`  | ``R^2 = \frac{\left( S_{XY} - \frac{1}{n} S_X S_Y \right)^2}{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}``  | True       | ``\mathcal L(x, y) = (x^2, xy, y^2, x, y)``               |
| Correlation Coefficient          | `:cor` | ``r = \frac{S_{XY} - \frac{1}{n} S_X S_Y}{\sqrt{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}}``             | True       | ``\mathcal L(x, y) = (x^2, xy, y^2, x, y)``               |
| Total Variation                  | `:tv`  | ``TV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) + \left( S_{YY} - \frac{1}{k} S_Y^2 \right)``                                                      | False      | ``\mathcal L(x, y) = (x^2, y^2, x, y)``                   |
| Covariance                       | `:cov` | ``\text{cov}(X, Y) = S_{XY} - \frac{1}{k} S_X S_Y``                                                                                                 | True       | ``\mathcal L(x, y) = (x, y, xy)``                         |
| Difference of Variances          | `:dv`  | ``DV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) - \left( S_{YY} - \frac{1}{k} S_Y^2 \right)``                                                      | True       | ``\mathcal L(x, y) = (x^2, y^2, x, y)``                   |
| Fraction of Variance Unexplained | `:fvu` | ``FVU = \frac{S_{YY}}{\left( S_{XX} - \frac{1}{k} S_X^2 \right)}``                                                                                  | False      | ``\mathcal L(x, y) = (x, x^2, y^2)``                      |

Note that some are minimized by default.
