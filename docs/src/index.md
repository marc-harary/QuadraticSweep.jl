# QuadraticSweep.jl Documentation

## Overview

**QuadraticSweep.jl** is a Julia package that provides efficient algorithms for selecting a subset of data points in 2D that maximizes various statistical criteria such as the coefficient of determination ($R^2$), correlation coefficient, total variation, etc. The package leverages combinatorial geometry, streaming algorithms, breadth-first search, and other optimization techniques to find these subsets.

The package implements two key algorithms:
- **Quadratic sweep algorithm**: Efficient topological sweep that projects data into higher dimensions to identify optimal subsets.
- **Brute-force combinatorial search**: Exhaustive search through all possible subsets.

## Installation

To install QuadraticSweep.jl, you can use Julia's package manager. This package is not yet registered in the General registry, so it needs to be installed directly from the GitHub repository:

```julia
using Pkg
Pkg.add(url="https://github.com/marc-harary/QuadraticSweep.jl")
```

## Main methods
```@docs
sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}
brute_force(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}
```

## Supported score functions
QuadraticSweep supports several score functions, and you can specify the score function to optimize by passing the appropriate symbol to the `sweep` function:

| Name                            | Symbol | Score Equation                                                                                                                                      | Maximizing | Embedding Equation                                       |
|----------------------------------|--------|-----------------------------------------------------------------------------------------------------------------------------------------------------|------------|----------------------------------------------------------|
| Coefficient of Determination     | `:r2`  | ``\frac{\left( S_{XY} - \frac{1}{n} S_X S_Y \right)^2}{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}``  | True       | ``(x^2, xy, y^2, x, y)``               |
| Correlation Coefficient          | `:cor` | ``\frac{S_{XY} - \frac{1}{n} S_X S_Y}{\sqrt{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}}``             | True       | ``(x^2, xy, y^2, x, y)``               |
| Total Variation                  | `:tv`  | ``\left( S_{XX} - \frac{1}{k} S_X^2 \right) + \left( S_{YY} - \frac{1}{k} S_Y^2 \right)``                                                      | False      | ``(x^2, y^2, x, y)``                   |
| Covariance                       | `:cov` | ``S_{XY} - \frac{1}{k} S_X S_Y``                                                                                                 | True       | ``\mathcal L(x, y) = (x, y, xy)``                         |
| Difference of Variances          | `:dv`  | ``\left( S_{XX} - \frac{1}{k} S_X^2 \right) - \left( S_{YY} - \frac{1}{k} S_Y^2 \right)``                                                      | True       | ``(x^2, y^2, x, y)``                   |
| Fraction of Variance Unexplained | `:fvu` | ``\frac{S_{YY}}{\left( S_{XX} - \frac{1}{k} S_X^2 \right)}``                                                                                  | False      | ``(x, x^2, y^2)``                      |

Note that some are minimized by default. 
