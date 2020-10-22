
struct Bounds
    lb
    ub
    
    function Bounds(lb, ub)
        lb = isnothing(lb) ? () : lb
        ub = isnothing(ub) ? () : ub

        size(lb) == () && size(ub) == () && !isnothing(lb) && !isnothing(ub) && return new(lb, ub)
        
        sizes = map(x -> size(x), [lb, ub])
        nonempty_bounds = findall(sizes .!= Ref(()))
        @assert !isempty(nonempty_bounds)
        length = sizes[nonempty_bounds[1]]
        
        lb = 1 in nonempty_bounds ? lb : repeat([-Inf], length)
        ub = 2 in nonempty_bounds ? ub : repeat([-Inf], length)

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
    δ
    state_bounds
    control_bounds
end

@with_kw mutable struct UninitializedProblem <: PSProblem
    dynamics
    path_constraints
    cost
    n_states
    n_controls
    n
    δ
    state_bounds::Bounds
    control_bounds::Bounds
    end_points
end



UninitializedProblem(prob::InitializedProblem) = UninitializedProblem(prob.dynamics, prob.path_constraints, prob.cost, size(prob.x₀)[1], size(prob.u₀)..., prob.δ, prob.end_points)

function symbolic_to_str(symbolic_expr::Basic; extend_dims=true, interpolate_variables=false)
    str_sym = string(symbolic_expr)
    
    # TODO: make it generic over variable names
    regex_dims = r"(\w+)(\d)\[(\d+)\]"
    regex_slice = r"([X|U]\[.*?\])"

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


function ode_constraint_exprs(M, dynamics, n, δ, ns, nc)

    function ode_error(dynamics, X, U)
        X * M' .- dynamics(X, U)
    end
    
    X_symbolic = zeros(Basic, (ns, n + 1))
    U_symbolic = zeros(Basic, (nc, n + 1))
    for (is, ic) in zip(1:ns, 1:nc), i in 1:n + 1
        X_symbolic[is, i] = symbols("X$is[$i]")[1]
        U_symbolic[ic, i] = symbols("U$ic[$i]")[1]
    end

    err = ode_error(dynamics, X_symbolic, U_symbolic)
    strs = symbolic_to_str.(err; interpolate_variables=false)
    strs .*= "≤ $δ"
    
    Meta.parse.(strs)
end

arr2d_to_arrlist(arr) = getindex.(Ref(arr), :, 1:size(arr)[2])
arrlist_to_arr2d(arr) = hcat(arr...)

function build_problem(prob::UninitializedProblem)
    n = prob.n
    τ, w = gausslobatto(n + 1)

    M = zeros(Float64, (n + 1, n + 1))
    diff_matrix!(M, n)

    xb = prob.state_bounds
    ub = prob.control_bounds
    broadcast_dynamics(X, U) = arrlist_to_arr2d(prob.dynamics.(arr2d_to_arrlist.([X, U])...))
    exprs = ode_constraint_exprs(M, broadcast_dynamics, n, prob.δ, prob.n_states, prob.n_controls)

    model = Model(Ipopt.Optimizer)
    
    @variables model begin
        xb.lb[i] <= X[i=1:prob.n_states, 1:n + 1] <= xb.ub[i], (base_name = "X$i")
        ub.lb[i] <= U[i=1:prob.n_controls, 1:n + 1] <= ub.ub[i], (base_name = "U$i")
    end

    @objective model Min prob.cost(X, U)
    @constraint model path_constraints 0 .<= prob.path_constraints(X, U)
    @constraint model initial_point X[:, 1] .== prob.end_points[1]
    @constraint model end_point X[:, end] .== prob.end_points[2]
    

    vars = Dict()
    for (is, ic) in zip(1:prob.n_states, 1:prob.n_controls), j in 1:n + 1
        vars[Meta.parse("X[$is, $j]")] = X[is, j]
        vars[Meta.parse("U[$ic, $j]")] = U[ic, j]
    end

    spliced_exprs = substitute_args.(exprs, Ref(vars))

    for expr in spliced_exprs
        add_NL_constraint(model, expr)
    end

    optimize!(model)
    
    X_sol = value.(X)
    U_sol = value.(U)

    cost_sol = prob.cost(X_sol, U_sol)

    X_sol, U_sol, cost_sol, model, τ
end