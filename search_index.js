var documenterSearchIndex = {"docs":
[{"location":"#QuadraticSweep.jl-Documentation","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"","category":"section"},{"location":"#Overview","page":"QuadraticSweep.jl Documentation","title":"Overview","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"QuadraticSweep.jl is a Julia package that provides efficient algorithms for selecting a subset of data points in 2D that maximizes various statistical criteria such as the coefficient of determination (R^2), correlation coefficient, total variation, etc. The package leverages combinatorial geometry, streaming algorithms, breadth-first search, and other optimization techniques to find these subsets.","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"The package implements two key algorithms:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Quadratic sweep algorithm: Efficient topological sweep that projects data into higher dimensions to identify optimal subsets.\nBrute-force combinatorial search: Exhaustive search through all possible subsets.","category":"page"},{"location":"#Installation","page":"QuadraticSweep.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"The package is not yet registered in the General registry, so it needs to be installed directly from the GitHub repository:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"using Pkg\nPkg.add(url=\"https://github.com/marc-harary/QuadraticSweep.jl\")","category":"page"},{"location":"#Main-methods","page":"QuadraticSweep.jl Documentation","title":"Main methods","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}\nbrute_force(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}","category":"page"},{"location":"#QuadraticSweep.sweep-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}}} where T<:Number","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.sweep","text":"sweep(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}\n\nEfficiently identifies the subset of k points from the dataset (x, y) that maximizes a given score function, using a topological sweep based on a lifting transformation.\n\nArguments\n\nx::Vector{T}: A vector of independent variable values, where T is a subtype of Number.\ny::Vector{T}: A vector of dependent variable values, with the same length as x.\nk::Int: The number of points to select in the optimal subset. Must be between 1 and n-1 where n is the total number of points.\nscore::Symbol: A symbol representing the score function to use (e.g., :r2, :cor, :tv). The symbol should correspond to a valid function in SCORE_FUNCTIONS.\n\nReturns\n\nA tuple containing:\nbest_idxs::Vector{Int64}: The indices of the k points that form the best subset.\nbest_score::Float64: The maximum score value for the best subset of size k.\n\nMethodology\n\nThe function applies a lifting transformation to the data points, projecting them into a higher-dimensional space. This facilitates the identification of an optimal separating hyperplane.\nIt then performs a topological sweep, iterating over candidate subsets (k-sets) to efficiently search for the subset with the highest score.\nIf the dataset size n is less than or equal to the embedding dimension d (from the score function), or if k <= d, the function falls back to a brute-force search.\n\nNotes\n\nThe score function and the corresponding lifting transformation are retrieved based on the provided score symbol.\nThe method ensures efficiency by skipping a brute-force search when possible, instead leveraging the topological sweep to reduce the number of computations.\n\nThrows\n\nAssertionError if x and y do not meet the required constraints (k < n, x and y are the same size).\nErrorException if the score symbol is invalid or not available in SCORE_FUNCTIONS.\n\nExample\n\ninlier_indices, inlier_cor = sweep(rand(25), rand(25); k = 20, score = :cor)\n\n\n\n\n\n","category":"method"},{"location":"#QuadraticSweep.brute_force-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}}} where T<:Number","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.brute_force","text":"brute_force(x::Vector{T}, y::Vector{T}; k::Int, score::Symbol) where {T <: Number}\n\nPerforms a brute-force search to find the subset of k points from the dataset (x, y) that maximizes the chosen score function.\n\nArguments\n\nx::Vector{T}: A vector of independent variable values, where T is a subtype of Number.\ny::Vector{T}: A vector of dependent variable values, with the same length as x.\nk::Int: The number of points to select in the optimal subset. Must be between 1 and n-1 where n is the total number of points.\nscore::Symbol: A symbol representing the score function to use (e.g., :r2, :cor, :tv). The symbol should correspond to a valid function in SCORE_FUNCTIONS.\n\nReturns\n\nA tuple containing:\nbest_idxs::Vector{Int64}: The indices of the k points that form the best subset.\nbest_score::Float64: The maximum score value for the best subset of size k.\n\nMethodology\n\nBrute-force combinatorial search.\n\nNotes\n\nExtremely slow for large n.\nUsed only as fallback method for sweep and for demonstration purposes.\n\nThrows\n\nAssertionError if x and y do not meet the expected constraints.\nErrorException if score is not valid.\n\nExample\n\ninlier_indices, inlier_tv = brute_force(rand(15), rand(15); k = 6, score = :tv)\n\n\n\n\n\n","category":"method"},{"location":"#Supported-score-functions","page":"QuadraticSweep.jl Documentation","title":"Supported score functions","text":"","category":"section"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"QuadraticSweep supports several score functions, and you can specify the score function to optimize by passing the appropriate symbol to the sweep function:","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Name Symbol Score Equation Maximizing Embedding Equation\nCoefficient of Determination :r2 fracleft( S_XY - frac1n S_X S_Y right)^2left( S_XX - frac1n S_X^2 right) left( S_YY - frac1n S_Y^2 right) True (x^2 xy y^2 x y)\nCorrelation Coefficient :cor fracS_XY - frac1n S_X S_Ysqrtleft( S_XX - frac1n S_X^2 right) left( S_YY - frac1n S_Y^2 right) True (x^2 xy y^2 x y)\nTotal Variation :tv left( S_XX - frac1k S_X^2 right) + left( S_YY - frac1k S_Y^2 right) False (x^2 y^2 x y)\nCovariance :cov S_XY - frac1k S_X S_Y True (x y xy)\nDifference of Variances :dv left( S_XX - frac1k S_X^2 right) - left( S_YY - frac1k S_Y^2 right) True (x^2 y^2 x y)\nFraction of Variance Unexplained :fvu fracS_YYleft( S_XX - frac1k S_X^2 right) False (x x^2 y^2)","category":"page"},{"location":"","page":"QuadraticSweep.jl Documentation","title":"QuadraticSweep.jl Documentation","text":"Note that some are minimized by default. ","category":"page"}]
}
