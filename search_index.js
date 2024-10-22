var documenterSearchIndex = {"docs":
[{"location":"#QuadraticSweep.jl-Documentation","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"","category":"section"},{"location":"#Overview","page":"QuadraticSweep.jl Documentation","title":"Overview","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"QuadraticSweep.jl is a Julia package that provides efficient algorithms for selecting a subset of data points in 2D that maximizes various statistical criteria such as the coefficient of determination (R^2), correlation coefficient, total variation, and more. The package leverages combinatorial geometry and optimization techniques to find these subsets.","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"The package implements two key algorithms:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Sweep Algorithm: Efficient topological sweep that projects data into higher dimensions to identify optimal subsets.\nBrute Force Algorithm: Exhaustive search through all possible subsets.","category":"page"},{"location":"#Installation","page":"QuadraticSweep.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"To install QuadraticSweep.jl, you can use Julia's package manager. This package is not yet registered in the General registry, so you'll need to install it directly from the GitHub repository:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"using Pkg\nPkg.add(url=\"https://github.com/username/QuadraticSweep.jl\")","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Ensure that you have Julia installed on your system. If not, download and install it from https://julialang.org/downloads/.","category":"page"},{"location":"#Usage","page":"QuadraticSweep.jl Documentation","title":"Usage","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Once installed, you can start using the QuadraticSweep module. Below is an example of how to use the sweep and brute_force functions to find the best subset of points maximizing a given score:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"using QuadraticSweep\n\n# Sample data\nx = [1.0, 2.0, 3.0, 4.0]\ny = [2.0, 3.0, 6.0, 8.0]\n\n# Finding the best subset using the R^2 score\nk = 2\nscore = :r2\nbest_idxs, best_score = sweep(x, y; k=k, score=score)\n\nprintln(\"Best subset indices: \", best_idxs)\nprintln(\"Best score: \", best_score)","category":"page"},{"location":"#Supported-score-functions","page":"QuadraticSweep.jl Documentation","title":"Supported score functions","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"QuadraticSweep supports several score functions, and you can specify the score function to maximize by passing the appropriate symbol to the sweep function: | Name                            | Symbol | Score Equation                                                                                                       | Maximizing | Embedding Equation                       | |––––––––––––––––-|––––|–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––|––––––|–––––––––––––––––––––| | Coefficient of Determination    | :r2  | R^2 = fracleft( S_XY - frac1n S_X S_Y right)^2left( S_XX - frac1n S_X^2 right) left( S_YY - frac1n S_Y^2 right) | True       | mathcal L(x y) = left( x^2 xy y^2 x y right) | | Correlation Coefficient         | :cor | r = fracS_XY - frac1n S_X S_Ysqrtleft( S_XX - frac1n S_X^2 right) left( S_YY - frac1n S_Y^2 right) | True       | mathcal L(x y) = left( x^2 xy y^2 x y right) | | Total Variation                 | :tv  | TV = left( S_XX - frac1k S_X^2 right) + left( S_YY - frac1k S_Y^2 right)                         | False      | mathcal L(x y) = left( x^2 y^2 x y right)     | | Covariance                      | :cov | textcov(X Y) = S_XY - frac1k S_X S_Y                                                                    | True       | mathcal L(x y) = left( x y xy right)           | | Difference of Variances         | :dv  | DV = left( S_XX - frac1k S_X^2 right) - left( S_YY - frac1k S_Y^2 right)                         | True       | mathcal L(x y) = left( x^2 y^2 x y right)     | | Fraction of Variance Unexplained| :fvu | FVU = fracS_YYleft( S_XX - frac1k S_X^2 right)                                                     | False      | mathcal L(x y) = left( x x^2 y^2 right)        |","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"func\nsweep","category":"page"},{"location":"#QuadraticSweep.func","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.func","text":"func(x)\n\nReturn double the number x plus 1.\n\n\n\n\n\n","category":"function"},{"location":"#QuadraticSweep.sweep","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.sweep","text":"sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}\n\nEfficiently identifies the subset of k points from the dataset (x, y) that maximizes a given score function, using a topological sweep based on a lifting transformation.\n\nArguments\n\nx::Vector{T}: A vector of independent variable values, where T is a subtype of Number.\ny::Vector{T}: A vector of dependent variable values, with the same length as x.\nk::Int: The number of points to select in the optimal subset. Must be between 1 and n-1 where n is the total number of points.\nscore::Symbol: A symbol representing the score function to use (e.g., :r2, :cor, :tv). The symbol should correspond to a valid function in SCORE_FUNCTIONS.\n\nReturns\n\nA tuple containing:\nbest_idxs::Vector{Int64}: The indices of the k points that form the best subset.\nbest_score::Float64: The maximum score value for the best subset of size k.\n\nMethodology\n\nThe function applies a lifting transformation to the data points, projecting them into a higher-dimensional space. This facilitates the identification of an optimal separating hyperplane.\nIt then performs a topological sweep, iterating over candidate subsets (k-sets) to efficiently search for the subset with the highest score.\nIf the dataset size n is less than or equal to the embedding dimension d (from the score function), or if k <= d, the function falls back to a brute-force search.\n\nNotes\n\nThe function uses a combination of combinatorics and linear algebra to determine the optimal subset.\nThe score function and the corresponding lifting transformation are retrieved based on the provided score symbol.\nThe method ensures efficiency by skipping a brute-force search when possible, instead leveraging the topological sweep to reduce the number of computations.\n\nThrows\n\nAssertionError if x and y do not meet the required constraints (k < n, x and y are the same size).\nErrorException if the score symbol is invalid or not available in SCORE_FUNCTIONS.\n\nExample\n\n```julia bestidxs, bestscore = sweep([1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 6.0, 8.0]; k = 2, score = :r2)\n\nSee also\n\nbrute_force: The fallback method for an exhaustive search.\ncheck_input: Validates inputs and retrieves the score function configuration.\n\n\n\n\n\n","category":"function"}]
}