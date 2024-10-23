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
inlier_indices, inlier_r2 = sweep(rand(20), rand(20); k = 10, score = :r2)
```

## Documentation
The documentation for the QuadraticSweep.jl package can be found [here](https://marc-harary.github.io/QuadraticSweep.jl/).
