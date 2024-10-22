# QuadraticSweep

[![arXiv](https://img.shields.io/badge/arXiv-2410.09316-b31b1b.svg)](https://arxiv.org/abs/2410.09316)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://marc-harary.github.io/QuadraticSweep.jl/)

[**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316)

The quadratic sweep algorithm deterministically detects optimal outliers for correlation analysis in efficient time.

## Documentation

The documentation for the QuadraticSweep.jl package can be found [here](https://marc-harary.github.io/QuadraticSweep.jl/).

## Usage
### Example

```julia
using QuadraticSweep
x, y = rand(10), rand(10)
inlier_idxs, inlier_score = sweep(x, y; k = 10, score = :r2)
```

### Score functions
| Name                            | Symbol | Score Equation                                                                                                       | Maximizing | Embedding Equation                       |
|---------------------------------|--------|----------------------------------------------------------------------------------------------------------------------|------------|------------------------------------------|
| Coefficient of Determination    | `:r2`  | $R^2 = \frac{\left( S_{XY} - \frac{1}{n} S_X S_Y \right)^2}{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}$ | True       | $\mathcal L(x, y) = \left( x^2, xy, y^2, x, y \right)$ |
| Correlation Coefficient         | `:cor` | $r = \frac{S_{XY} - \frac{1}{n} S_X S_Y}{\sqrt{\left( S_{XX} - \frac{1}{n} S_X^2 \right) \left( S_{YY} - \frac{1}{n} S_Y^2 \right)}}$ | True       | $\mathcal L(x, y) = \left( x^2, xy, y^2, x, y \right)$ |
| Total Variation                 | `:tv`  | $TV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) + \left( S_{YY} - \frac{1}{k} S_Y^2 \right)$                         | False      | $\mathcal L(x, y) = \left( x^2, y^2, x, y \right)$     |
| Covariance                      | `:cov` | $\text{cov}(X, Y) = S_{XY} - \frac{1}{k} S_X S_Y$                                                                    | True       | $\mathcal L(x, y) = \left( x, y, xy \right)$           |
| Difference of Variances         | `:dv`  | $DV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) - \left( S_{YY} - \frac{1}{k} S_Y^2 \right)$                         | True       | $\mathcal L(x, y) = \left( x^2, y^2, x, y \right)$     |
| Fraction of Variance Unexplained| `:fvu` | $FVU = \frac{S_{YY}}{\left( S_{XX} - \frac{1}{k} S_X^2 \right)}$                                                     | False      | $\mathcal L(x, y) = \left( x, x^2, y^2 \right)$        |

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
