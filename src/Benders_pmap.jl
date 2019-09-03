mutable struct Solution
    feasible::Bool
    objective_value::Union{Nothing, Float64}
    # Map between the variable ref
    # and the primal value
    variable_value::Dict{JuMP.VariableRef, Float64}
    # Map between the parameter
    # and the dual value
    parameter_dual::Dict{ParameterRef, Float64}
end

function optimize(model::ParametrizedModel)
    JuMP.optimize!(model.model)
    feasible = JuMP.primal_status(model.model) == MOI.FEASIBLE_POINT
    if feasible
        @assert JuMP.termination_status(model.model) == MOI.OPTIMAL
        objective_value = JuMP.objective_value(model.model)
    else
        @assert JuMP.termination_status(model.model) in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
        if isempty(model.parameter_map)
            objective_value = nothing
        else
            if JuMP.termination_status(model.model) == MOI.ALMOST_INFEASIBLE
                @warn("Using infeasibility ray of result with termination status: ALMOST_INFEASIBLE.")
            end
            @assert JuMP.dual_status(model.model) == MOI.INFEASIBILITY_CERTIFICATE
            objective_value = JuMP.dual_objective_value(model.model)
        end
    end
    variable_value = Dict{JuMP.VariableRef, Float64}()
    if feasible
        for vref in values(model.variable_map)
            variable_value[vref] = JuMP.value(vref)
        end
        for θ in values(model.θ)
            variable_value[θ] = JuMP.value(θ)
        end
    end
    parameter_dual = Dict{ParameterRef, Float64}()
    for parameter in values(model.parameter_map)
        parameter_dual[parameter] = JuMP.dual(parameter)
    end
    Solution(feasible, objective_value, variable_value, parameter_dual)
end

function set_parent_solution!(model::ParametrizedModel, parent::ParametrizedModel, parent_solution::Solution)
    for (index, parameter) in model.parameter_map
        vref = parent.variable_map[index]
        value = parent_solution.variable_value[vref]
        JuMP.fix(parameter, value)
    end
end

function add_cutting_planes(master_model, master_solution, sub_models, sub_solutions, TOL)
    cut_added = false
    # add cutting planes, one per scenario
    for (id, sol) = sub_solutions
        sub_model = sub_models[id]
        aff = JuMP.AffExpr(sol.objective_value)
        for (index, parameter) in sub_model.parameter_map
            dual = sol.parameter_dual[parameter]
            vref = master_model.variable_map[index]
            aff.constant -= dual * master_solution.variable_value[vref]
            JuMP.add_to_expression!(aff, dual, vref)
        end
        if sol.feasible
            JuMP.add_to_expression!(aff, -1.0, master_model.θ[id])
            # Check if the cut is useful
            if JuMP.value(aff, vref -> master_solution.variable_value[vref]) - TOL < 0
                continue
            end
        else
            # It is a ray, so `alpha * aff` is also valid for any alpha >= 0.
            # Hence `aff` might have very large coefficients and alter
            # the numerical accuracy of the master's solver.
            # We scale it to avoid this issue
            scaling = abs(aff.constant)
            if iszero(scaling)
                scaling = maximum(term -> term[1], JuMP.linear_terms(aff))
            end
            scaled_aff = typeof(aff)(sign(aff.constant))
            for (coef, var) in JuMP.linear_terms(aff)
                JuMP.add_to_expression!(scaled_aff, coef / scaling, var)
            end
            aff = scaled_aff
        end
        @constraint(master_model.model, aff <= 0)
        cut_added = true
    end
    return cut_added
end

function Benders_pmap(master_model::ParametrizedModel,
                      sub_models::Dict{Int, ParametrizedModel},
                      TOL=1e-5)::Solution
    master_solution = nothing
    cut_added = true
    while cut_added
        master_solution = optimize(master_model)
        if !master_solution.feasible
            break
        end

        for (id, sub_model) in sub_models
            set_parent_solution!(sub_model, master_model, master_solution)
        end

        sub_solutions = Dict(map(x -> (x[1] => optimize(x[2])), collect(sub_models))...)
        cut_added = add_cutting_planes(master_model, master_solution, sub_models, sub_solutions, TOL)
    end
    return master_solution
end
