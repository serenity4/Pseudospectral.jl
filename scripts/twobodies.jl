using JuMP
using ApproxFun
using Pseudospectral
using Plots
using FastGaussQuadrature


const M = zeros(n + 1, n + 1) ; diff_matrix!(M, n)
const w = gausslobatto(n + 1)[2]

function cost(X, U)
    integrand = w .* ((M * U[2, :]).^2)# .+ 0.01 .* (M * U[1, :]).^2)
    sum(integrand)
end

function dynamics(X, U)
    [cos(U[1]) * U[2], sin(U[1]) * U[2]]
end

function path_constraints(X, U)
    sum((X .- 1).^2, dims=1) .- 0.6^2
end

function circle_shape(radius; origin=(0., 0.), n=50)
    xs = range(0, 2pi * (1 + 1 / n), step=2pi / n) |> x -> radius * sin.(x) |> Array{Float64}
    ys = range(0, 2pi * (1 + 1 / n), step=2pi / n) |> x -> radius * cos.(x) |> Array{Float64}
    origin[1] .+ xs, origin[2] .+ ys
end

n = 50

# TODO: allow for non-uniform bounds
prob = UninitializedProblem(
    dynamics,
    path_constraints,
    cost,
    2,
    2,
    n,
    1e-6,
    Bounds([0., 0.], [2., 2.]),
    Bounds([0, 0], [2pi, 2.]),
    ([0., 0.], [2., 2.])
)


X, U, final_cost, model, τ = Pseudospectral.build_problem(prob);
p = plot(X[1, :], X[2, :], marker_z=τ, marker=:o)
plot!(p, circle_shape(0.6, origin=(1., 1.))..., label="obstacle")
savefig(p, "res.png")
scatter(τ, U[1, :], marker_z=τ, title="θ")
scatter(τ, U[2, :], marker_z=τ, title="v")

termination_status(model)
# primal_status(model)
# dual_status(model)


using DataInterpolations


function reconstruct(X, τ; step=0.01)
    map.(LagrangeInterpolation.([X[1, :], X[2, :], τ], Ref(τ), Ref(length(τ) - 1)), Ref(-1:step:1))
end

lag_x, lag_y = LagrangeInterpolation.([X[1, :], X[2, :]], Ref(τ), Ref(n))
lag_τ = LagrangeInterpolation(τ, τ)


plot(X[1, :], X[2, :], marker_z=lag_τ, marker=:o)
scatter(X[1, :], X[2, :], marker_z=τ)
τ_subsampled = -1:0.01:1

xs, ys, ts = reconstruct(X, τ)
scatter(xs, ys, marker_z=ts)

