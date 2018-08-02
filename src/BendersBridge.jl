using ParameterJuMP
using MathOptInterface
const MOI = MathOptInterface

export BendersBridge

function parametrized_function(var::StructuredVariableRef,
                               structured_model::StructuredModel,
                               model::JuMP.Model,
                               variable_map::Dict{Int, JuMP.VariableRef},
                               parameter_map::Dict{Int, ParameterJuMP.Parameter})
    if var.model === structured_model
        return variable_map[var.idx]
    elseif var.model === structured_model.parent
        if !haskey(parameter_map, var.idx)
            parameter_map[var.idx] = ParameterJuMP.Parameter(model)
        end
        return parameter_map[var.idx]
    else
        error("A StructuredModel cannot have expression using variables of a StructuredModel different from itself or its parent")
    end
end
function parametrized_function(aff::JuMP.GenericAffExpr{C, StructuredVariableRef},
                               structured_model::StructuredModel,
                               model::JuMP.Model,
                               variable_map::Dict{Int, JuMP.VariableRef},
                               parameter_map::Dict{Int, ParameterJuMP.Parameter}) where C
    param_aff = ParameterJuMP.PAE{C}(JuMP.GenericAffExpr{C, JuMP.VariableRef}(aff.constant),
                                     JuMP.GenericAffExpr{C, ParameterJuMP.Parameter}(zero(C)))
    for (coef, var) in JuMP.linearterms(aff)
        JuMP.add_to_expression!(param_aff,
                                coef,
                                parametrized_function(var, structured_model, model, variable_map, parameter_map))
    end
    # The inferred return value is an union of two concrete types but
    # small unions are efficiently handled in Julia v0.7
    if JuMP.iszero(param_aff.p)
        return param_aff.v
    else
        return param_aff
    end
end
function parametrized_function(quad::JuMP.GenericQuadExpr{C, StructuredVariableRef},
                               structured_model::StructuredModel,
                               model::JuMP.Model,
                               variable_map::Dict{Int, JuMP.VariableRef},
                               parameter_map::Dict{Int, ParameterJuMP.Parameter}) where C
    param_aff = parametrized_function(quad.aff, structured_model, model, variable_map, parameter_map)
    if param_aff isa ParameterJuMP.PAE
        error("parametrized quadratic functions are not supported yet")
        #quadv = JuMP.GenericQuadExpr{C, JuMP.VariableRef}(param_aff.v)
        #quadp = JuMP.GenericQuadExpr{C, ParameterJuMP.Parameter}(param_aff.p)
    else
        quadv = JuMP.GenericQuadExpr{C, JuMP.VariableRef}(param_aff)
        quadp = JuMP.GenericQuadExpr{C, ParameterJuMP.Parameter}(zero(C))
    end
    for (coef, var1, var2) in JuMP.quadterms(quad)
        if var1.model === structured_model && var2.model == structured_model
            JuMP.add_to_expression!(quadv, coef, variable_map[var1.idx], variable_map[var2.idx])
        else
            error("parametrized quadratic functions are not supported yet")
        end
    end
    if JuMP.iszero(quadp)
        return quadv
    else
        error("parametrized quadratic functions are not supported yet")
    end
end
parametrized_function(funcs::Vector, args...) = map(f -> parametrized_function(f, args...), funcs)

struct ParametrizedModel
    structured_model::StructuredModel
    model::JuMP.Model
    # Map between index of structured variable in `structured_model` and
    # the corresponding variable in `model`.
    variable_map::Dict{Int, JuMP.VariableRef}
    # Map between index of structured variable in `structured_model.parent`
    # and the corresponding parameter in `model`.
    parameter_map::Dict{Int, ParameterJuMP.Parameter}
    # Cost of children
    θ::Dict{Int, JuMP.VariableRef}
end
function ParametrizedModel(structured_model::StructuredModel, args...; kwargs...)
    if structured_model.parent === nothing
        model = Model(args...; kwargs...)
    else
        model = ModelWithParams(args...; kwargs...)
    end
    variable_map = Dict{Int, JuMP.VariableRef}()
    for (index, var) in structured_model.variables
        name = structured_model.varnames[index]
        variable_map[index] = JuMP.addvariable(model, var, name)
    end
    parameter_map = Dict{Int, ParameterJuMP.Parameter}()
    for (index, con) in structured_model.constraints
        name = structured_model.connames[index]
        param_fun = parametrized_function(con.func, structured_model, model,
                                          variable_map, parameter_map)
        param_con = JuMP.buildconstraint(error, param_fun, con.set)
        JuMP.addconstraint(model, param_con, name)
    end
    objective_function = parametrized_function(structured_model.objective_function,
                                               structured_model, model,
                                               variable_map, parameter_map)
    if structured_model.objective_sense == :Max
        objective_function = -objective_function
    end
    θ = Dict{Int, VariableRef}()
    for (id, proba) in getprobability(structured_model)
        θid = @variable(model, lowerbound = 0.0)
        θ[id] = θid
        JuMP.add_to_expression!(objective_function, proba, θid)
    end
    JuMP.setobjective(model, structured_model.objective_sense, objective_function)
    ParametrizedModel(structured_model, model, variable_map, parameter_map, θ)
end

include("Benders_pmap.jl")

# The optimizers are function returning an empty optimizer,
# this will be replaced by factories once this is implemented in JuMP
function BendersBridge(structured_model::StructuredModel, master_optimizer::Function, sub_optimizer::Function)
    master_model = ParametrizedModel(structured_model, optimizer=master_optimizer())
    children = getchildren(structured_model)
    sub_models = Dict{Int, ParametrizedModel}()
    for (id, child) in getchildren(structured_model)
        sub_models[id] = ParametrizedModel(child, optimizer=sub_optimizer())
    end
    return Benders_pmap(master_model, sub_models)
end
