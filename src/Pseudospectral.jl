module Pseudospectral

using FastGaussQuadrature, ApproxFun, Parameters, JuMP, SymEngine, Ipopt

include("diff_matrix.jl")
include("problem.jl")

export diff_matrix!, Bounds, UninitializedProblem, PSProblem, build_problem

end # module
