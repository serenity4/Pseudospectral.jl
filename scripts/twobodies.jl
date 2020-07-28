using JuMP
using ApproxFun
using Pseudospectral
using Plots
using FastGaussQuadrature



function dynamics(X, U)
    res = @. [sin(U[1, :]) * U[2, :], cos(U[1, :]) * U[2, :]]
    transpose(hcat(res...))
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
prob = UninitializedProblem(
    dynamics,
    path_constraints,
    cost,
    2,
    2,
    n,
    1e-6,
    # Bounds([0., 0.], [2., 2.]),
    # Bounds([0., 0.], [2pi, 2pi]),
    Bounds(0., 2.),
    Bounds(0, 2.),
    ([0., 0.], [2., 2.])
)

const M = zeros(n + 1, n + 1) ; diff_matrix!(M, n)
const w = gausslobatto(n + 1)[2]
function cost(X, U)
    sum(w .* (M * U[1]).^2 + 0.1 * (M * (U[1] .- pi / 4).^2))
end

X, U, final_cost, model, τ = build_problem(prob);
p = plot(X[1, :], X[2, :], marker_z=τ, marker=:o)
plot!(p, circle_shape(0.6, origin=(1., 1.))..., label="obstacle")
savefig(p, "res.png")
scatter(τ, θ, marker_z=τ)
scatter(τ, v, marker_z=τ)

termination_status(model)
# primal_status(model)
# dual_status(model)


using DataInterpolations


lag_x, lag_y = LagrangeInterpolation.([x[1, :], x[2, :]], Ref(τ), Ref(n))
lag_τ = LagrangeInterpolation(τ, τ)


plot(x[1, :], x[2, :], marker_z=lag_τ, marker=:o)
scatter(x[1, :], x[2, :], marker_z=τ)
# scatter(lag_x, lag_y, marker_z = lag_τ)

