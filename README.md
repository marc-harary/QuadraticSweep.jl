# QuadraticSweep

[![arXiv](https://img.shields.io/badge/arXiv-2410.09316-b31b1b.svg)](https://arxiv.org/abs/2410.09316)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://marc-harary.github.io/QuadraticSweep.jl/)

Official implementation of the quadratic sweep algorithm from [**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316). Subsets of data that maximize statistical criteria like $R^2$ and correlation coefficient are optimally identified in efficient time via a topological sweep in higher-dimensional embedding spaces.

## Installation

Currently, the package is not added to the general Julia registry; it will instead need to be installed directly from this repository:

```julia
using Pkg
Pkg.add(url="https://github.com/marc-harary/QuadraticSweep.jl")
```

## Usage

Calling the core `sweep` method is straightforward. Simply pass in the `x` and `y` data vectors, the number `k` of total inliers, and the desired `score` function:

```julia
using QuadraticSweep
x, y = rand(20), rand(20)
inlier_indices, inlier_r2 = sweep(x, y; k = 17, score = :r2)
```

## Documentation

The documentation for the QuadraticSweep.jl package can be found [here](https://marc-harary.github.io/QuadraticSweep.jl/). It includes supported score functions and more detailed explanations of this package's utilities and core functionality.

## Citation

If you use this code, please cite:

```bibtex
@article{harary2024combinatorial,
  title={Combinatorial optimization of the coefficient of determination},
  author={Harary, Marc},
  journal={arXiv preprint arXiv:2410.09316},
  year={2024}
}
```
