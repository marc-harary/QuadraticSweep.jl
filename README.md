# QuadraticSweep

[![arXiv](https://img.shields.io/badge/arXiv-2410.09316-b31b1b.svg)](https://arxiv.org/abs/2410.09316)

[**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316)

The quadratic sweep algorithm deterministically detects optimal outliers for correlation analysis in efficient time.

## Usage
### Example

```julia
using QuadraticSweep
x, y = rand(10), rand(10)
inlier_idxs, inlier_score = sweep(x, y; k = 10, score = :r2)
```

### Score functions
| Name                          | Symbol | Score Equation                                                                                                    | Maximizing | Embedding Dimension | Embedding Equation                |
|-------------------------------|--------|-------------------------------------------------------------------------------------------------------------------|-----------------------|---------------------|-----------------------------------|
| Coefficient of Determination   | `:r2`  | \( R^2 = \frac{(S_{XY} - \frac{1}{n} S_X S_Y)^2}{(S_{XX} - \frac{1}{n} S_X^2) (S_{YY} - \frac{1}{n} S_Y^2)} \)    | True            | 5                   | \( (x^2, xy, y^2, x, y) \)       |
| Correlation Coefficient        | `:cor` | \( r = \frac{S_{XY} - \frac{1}{n} S_X S_Y}{\sqrt{(S_{XX} - \frac{1}{n} S_X^2)(S_{YY} - \frac{1}{n} S_Y^2)}} \)    | True            | 5                   | \( (x^2, xy, y^2, x, y) \)       |
| Total Variation                | `:tv`  | \( TV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) + \left( S_{YY} - \frac{1}{k} S_Y^2 \right) \)                   | False            | 4                   | \( (x^2, y^2, x, y) \)           |
| Covariance                     | `:cov` | \( Cov(X, Y) = S_{XY} - \frac{1}{k} S_X S_Y \)                                                                    | True            | 3                   | \( (x, y, xy) \)                 |
| Difference of Variances        | `:dv`  | \( DV = \left( S_{XX} - \frac{1}{k} S_X^2 \right) - \left( S_{YY} - \frac{1}{k} S_Y^2 \right) \)                   | True            | 4                   | \( (x^2, y^2, x, y) \)           |
| Fraction of Variance Unexplained| `:fvu` | \( FVU = \frac{(S_{YY}}{S_{XX} - \frac{1}{k} S_X^2} \)                           | False            | 3                   | \( (x, x^2, y^2) \)                 |

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
