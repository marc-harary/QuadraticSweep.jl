# QuadraticSweep

[![arXiv](https://img.shields.io/badge/arXiv-2410.09316-b31b1b.svg)](https://arxiv.org/abs/2410.09316)

[**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316)

The quadratic sweep algorithm deterministically detects optimal outliers for correlation analysis in efficient time.

## Build
`QuadraticSweep.jl` can be installed via
```julia
using Pkg
Pkg.add(url="https://github.com/marc-harary/QuadraticSweep.jl")
```
To only build the package locally, run
```bash
julia --project --color=yes -e 'using Pkg; Pkg.instantiate()'
```

## Tests/reproducibility
The tests included in the original manuscript can be performed by running
```julia
using Pkg
Pkg.test("QuadraticSweep")
```
if `QuadraticSweep.jl` is installed. If it is instead locally built, run
```bash
julia --project --color=yes -e 'using Pkg; Pkg.test()'
```
