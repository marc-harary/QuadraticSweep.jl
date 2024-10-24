# Experiments
The following notebooks reproduce the key experiments in [**Combinatorial optimization of the coefficient of determination**](https://arxiv.org/abs/2410.09316).
## Notebooks
The notebooks are as follows:
- `seeded_tests.ipynb`: Thoroughly tests the quadratic sweep against brute-force combinatorial search for all scoring functions.
- `separability_tests.ipynb`: Uses [Ipopt](https://github.com/jump-dev/Ipopt.jl) and [JuMP](https://github.com/jump-dev/JuMP.jl) to determine the quadratic separability of inlying and outlying subsets following different lifting functions.
- `heuristics_tests.ipynb`: Tests the performance of the quadratic sweep against various heuristic methods.
## Setup
You will need an environment specific to the `notebooks` directory with its specific dependencies. Run the following to interact with the notebooks:
```julia
using IJulia
# Install the Jupyter kernel specifically for this environment:
IJulia.installkernel("QuadraticSweep", "--project=$(Base.current_project())")
```
For `heuristics_tests.ipynb`, you will need to have Python installed along with Scikit-Learn and [simanneal](https://github.com/perrygeo/simanneal).
