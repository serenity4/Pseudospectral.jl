using JuMP
using Ipopt
using ApproxFun
using Pseudospectral
using FastGaussQuadrature
using SymEngine
using Plots

function symbolic_to_str(symbolic_expr::Basic; extend_dims=true, interpolate_variables=false)
    str_sym = string(symbolic_expr)
    
    # TODO: make it generic over variable names
    regex_dims = r"(\w+)(\d)\[(\d+)\]"
    regex_slice = r"([x|v|θ]\[.*?\])"

    if extend_dims
        str_sym = replace(str_sym, regex_dims => s"\1[\2,\3]")
    end
    if interpolate_variables
        str_sym = replace(str_sym, regex_slice => s"$(\1)")
    end

    str_sym
end

substitute_args(expr, vars) = expr
substitute_args(expr::Symbol, vars) = get(vars, expr, expr)
function substitute_args(expr::Expr, vars)
    for (i, arg) in enumerate(expr.args)
        if arg in keys(vars)
            expr.args[i] = vars[arg]
        else
            expr.args[i] = substitute_args(arg, vars)
        end
    end
    return expr
end


function ode_constraint_exprs(M, dynamics, n, δ)

    function ode_error(dynamics, x, θ, v)
        x * M' - hcat(dynamics.(θ, v)...)
    end
    
    x_symbolic = zeros(Basic, (2, n + 1))
    θ_symbolic = [symbols("θ[$i]") for i in 1:n + 1]
    v_symbolic = [symbols("v[$i]") for i in 1:n + 1]
    for i in 1:n + 1
        x_symbolic[1, i] = symbols("x1[$i]")[1]
        x_symbolic[2, i] = symbols("x2[$i]")[1]
    end

    err = ode_error(dynamics, x_symbolic, θ_symbolic, v_symbolic)
    strs = symbolic_to_str.(err; interpolate_variables=false)
    strs .*= "≤ $δ"
    
    Meta.parse.(strs)
end

function build_problem(n, δ)
    τ, w = gausslobatto(n + 1)

    M = zeros(Float64, (n + 1, n + 1))
    diff_matrix!(M, n)

    # function dynamics(X, U)
    #     [sin(U[1]) * U[2], cos(U[1]) * U[2]]
    # end

    # function cost(X, U)
    #     sum(w .* (M * U[1]).^2)
    # end

    function dynamics(θ, v)
        [sin(θ) * v, cos(θ) * v]
    end

    function cost(x, θ, v)
        sum(w .* ((M * v).^2 .+ 0.1 * (M * (θ .- pi / 4).^2)))
    end

    function path_constraints(x, θ, v)
        sum((x .- 1).^2, dims=1) .- 0.6^2
    end

    x0 = [0, 0]
    xf = [2, 2]
    
    xb = ([0, 0], [2, 2])
    θb = [-pi, pi]
    vb = [0, 2]
    
    model = Model(Ipopt.Optimizer)
    
    @variable model 0 <= x[1:2, 1:n + 1] <= 2 base_name = "position"
    # @variable model x[1:2, 1:n + 1] base_name = "position"
    
    @variable model 0 <= θ[1:n + 1] <= 2pi base_name = "rotation angle"
    # @variable model 0 <= v[1:n + 1] base_name = "speed"
    @variable model vb[1] ≤ v[1:n + 1] ≤ vb[2] base_name = "speed"
    
    @objective model Min cost(x, θ, v)
    
    @constraint model collision 0 .<= path_constraints(x, θ, v)
    @constraint model initial_point x[:, 1] .== x0
    @constraint model end_point x[:, end] .== xf

    
    exprs = ode_constraint_exprs(M, dynamics, n, δ)

    vars = Dict()
    for i in 1:2, j in 1:n + 1
        vars[Meta.parse("x[$i, $j]")] = x[i, j]
        vars[Meta.parse("v[$j]")] = v[j]
        vars[Meta.parse("θ[$j]")] = θ[j]
    end


    spliced_exprs = substitute_args.(exprs, Ref(vars))

    for expr in spliced_exprs
        add_NL_constraint(model, expr)
    end

    optimize!(model)
    
    x_sol = value.(x)
    θ_sol = value.(θ)
    v_sol = value.(v)

    cost_sol = cost(x_sol, θ_sol, v_sol)

    x_sol, θ_sol, v_sol, cost_sol, model, τ
end

function circle_shape(radius; origin=(0., 0.), n=50)
    xs = range(0, 2pi * (1 + 1 / n), step=2pi / n) |> x -> radius * sin.(x) |> Array{Float64}
    ys = range(0, 2pi * (1 + 1 / n), step=2pi / n) |> x -> radius * cos.(x) |> Array{Float64}
    origin[1] .+ xs, origin[2] .+ ys
end

n = 50
x, θ, v, cost, model, τ = build_problem(n, 1e-6);
p = plot(x[1, :], x[2, :], marker_z=τ, marker=:o)
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

