
struct Bounds
    lb
    ub
    
    function Bounds(lb, ub)
        lengths = map(x -> first(size(x)), lb, ub)
        lb = isnothing(lb) ? repeat([-Inf], lengths[1]) : lb
        ub = isnothing(lb) ? repeat([Inf], lengths[2]) : ub
        new(lb, ub)
    end
end

abstract type PSProblem end

@with_kw mutable struct InitializedProblem{T} <: PSProblem
    dynamics
    path_constraints
    cost
    x₀::AbstractArray{T}
    u₀::AbstractArray{T}
    δ = 1e-6
    state_bounds = Bounds(nothing, nothing)
    control_bounds = Bounds(nothing, nothing)
end

@with_kw mutable struct UnitializedProblem{T} <: PSProblem
    dynamics
    path_constraints
    cost
    n_states
    n_controls
    n
    δ
    state_bounds::Bounds
    control_bounds::Bounds
    end_points::Bounds
end



UnitializedProblem(prob::InitializedProblem) = UninitializedProblem(prob.dynamics, prob.path_constraints, prob.cost, size(prob.x₀)[1], size(prob.u₀)..., prob.δ, prob.end_points)


# function build_problem(prob::UninitializedProblem)
#     n = prob.n
#     τ, w = gausslobatto(n + 1)

#     M = zeros(Float64, (n + 1, n + 1))
#     diff_matrix!(M, n)

#     x0, xf = prob.end_points
#     xb = prob.state_bounds
#     ub = prob.control_bounds


#     θb = [-pi, pi]
#     vb = [0, 2]
    
#     model = Model(Ipopt.Optimizer)
    
#     @variable model 0 <= x[1:2, 1:n + 1] <= 2 base_name = "position"
#     # @variable model x[1:2, 1:n + 1] base_name = "position"
    
#     @variable model 0 <= θ[1:n + 1] <= 2pi base_name = "rotation angle"
#     # @variable model 0 <= v[1:n + 1] base_name = "speed"
#     @variable model vb[1] ≤ v[1:n + 1] ≤ vb[2] base_name = "speed"
    
#     @objective model Min cost(x, θ, v)
    
#     @constraint model collision 0 .<= path_constraints(x, θ, v)
#     @constraint model initial_point x[:, 1] .== x0
#     @constraint model end_point x[:, end] .== xf

    
#     exprs = ode_constraint_exprs(M, dynamics, n, δ)

#     vars = Dict()
#     for i in 1:2, j in 1:n + 1
#         vars[Meta.parse("x[$i, $j]")] = x[i, j]
#         vars[Meta.parse("v[$j]")] = v[j]
#         vars[Meta.parse("θ[$j]")] = θ[j]
#     end


#     spliced_exprs = substitute_args.(exprs, Ref(vars))

#     for expr in spliced_exprs
#         add_NL_constraint(model, expr)
#     end

#     optimize!(model)
    
#     x_sol = value.(x)
#     θ_sol = value.(θ)
#     v_sol = value.(v)

#     cost_sol = cost(x_sol, θ_sol, v_sol)

#     x_sol, θ_sol, v_sol, cost_sol, model, τ
# end