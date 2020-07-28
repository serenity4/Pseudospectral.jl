module Pseudospectral

using FastGaussQuadrature, ApproxFun, Parameters

include("diff_matrix.jl")
include("problem.jl")

export diff_matrix!, Bounds, UnitializedProblem, PSProblem

end # module
