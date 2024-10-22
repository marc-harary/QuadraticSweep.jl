# QuadraticSweep

[![arXiv](https://img.shields.io/badge/arXiv-2410.09316-b31b1b.svg)](https://arxiv.org/abs/2410.09316)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://marc-harary.github.io/QuadraticSweep.jl/)

Official implementation of the quadratic sweep algorithm from [**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316).

## Installation
```julia
using Pkg
Pkg.add(url="https://github.com/marc-harary/QuadraticSweep.jl")
```

## Usage
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

## Documentation
The documentation for the QuadraticSweep.jl package can be found [here](https://marc-harary.github.io/QuadraticSweep.jl/).
