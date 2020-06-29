module StructJuMP

export StructuredModel

using MathOptInterface
const MOI = MathOptInterface

using JuMP # To reexport, should be using (not import)
# Macro to exportall
macro exportall(pkg)
    Expr(:export, names(JuMP)...)
end
@exportall JuMP


# The following is largely inspired from JuMP/test/JuMPExtension.jl

mutable struct StructuredModel <: JuMP.AbstractModel
    # Structured data
    parent::Union{StructuredModel, Nothing}
    children::Dict{Int, StructuredModel}
    probability::Dict{Int, Float64}
    num_scen::Int

    # Model data
    nextvaridx::Int                                 # Next variable index is nextvaridx+1
    variables::Dict{Int, JuMP.AbstractVariable}     # Map varidx -> variable
    varnames::Dict{Int, String}                     # Map varidx -> name
    nextconidx::Int                                 # Next constraint index is nextconidx+1
    constraints::Dict{Int, JuMP.AbstractConstraint} # Map conidx -> variable
    connames::Dict{Int, String}                     # Map varidx -> name
    objective_sense::MOI.OptimizationSense
    objective_function::JuMP.AbstractJuMPScalar
    objdict::Dict{Symbol, Any}                      # Same that JuMP.Model's field `objdict`
    function StructuredModel(; parent=nothing, same_children_as=nothing, id=0,
                             num_scenarios::Int=0,
                             prob::Float64=default_probability(parent))
        if same_children_as !== nothing
            if !isa(same_children_as, StructuredModel)
                error("The JuMP model given for the argument `same_children_as' is not valid. Please create it using the `StructuredModel' function.")
            end
            probability = same_children_as.probability
            children = same_children_as.children
        else
            probability = Dict{Int, Float64}()
            children = Dict{Int, JuMP.Model}()
        end

        model = new(parent, children, probability, num_scenarios,                    # Structured
                    0, Dict{Int, JuMP.AbstractVariable}(),   Dict{Int, String}(),    # Model Variables
                    0, Dict{Int, JuMP.AbstractConstraint}(), Dict{Int, String}(),    # Model Constraints
                    MOI.FEASIBILITY_SENSE, zero(JuMP.GenericAffExpr{Float64, StructuredVariableRef}), # Model objective
                    Dict{Symbol, Any}())                                             # Model objects

        if parent === nothing
            id = 0
        else
            @assert id != 0
            parent.children[id] = model
            parent.probability[id] = prob
        end

        return model
    end
end

#### Structured ####

getparent(model::StructuredModel)      = model.parent
getchildren(model::StructuredModel)    = model.children
getprobability(model::StructuredModel) = model.probability
num_scenarios(model::StructuredModel)  = model.num_scen

default_probability(model::StructuredModel) = 1 / num_scenarios(model)
default_probability(::Nothing) = 1.0

#### Model ####

JuMP.object_dictionary(m::StructuredModel) = m.objdict

# Variables
struct StructuredVariableRef <: JuMP.AbstractVariableRef
    model::StructuredModel # `model` owning the variable
    idx::Int       # Index in `model.variables`
end
Base.broadcastable(v::StructuredVariableRef) = Ref(v)
Base.copy(v::StructuredVariableRef) = v
Base.:(==)(v::StructuredVariableRef, w::StructuredVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::StructuredVariableRef) = v.model
JuMP.isequal_canonical(v::StructuredVariableRef, w::StructuredVariableRef) = v == w
JuMP.variable_type(::StructuredModel) = StructuredVariableRef
function JuMP.add_variable(m::StructuredModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = StructuredVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.set_name(vref, name)
    vref
end
function MOI.delete!(m::StructuredModel, vref::StructuredVariableRef)
    delete!(m.variables, vref.idx)
    delete!(m.varnames, vref.idx)
end
MOI.is_valid(m::StructuredModel, vref::StructuredVariableRef) = vref.idx in keys(m.variables)
JuMP.num_variables(m::StructuredModel) = length(m.variables)

# Internal function
variable_info(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info
function update_variable_info(vref::StructuredVariableRef, info::JuMP.VariableInfo)
    vref.model.variables[vref.idx] = JuMP.ScalarVariable(info)
end

JuMP.has_lower_bound(vref::StructuredVariableRef) = variable_info(vref).has_lb
function JuMP.lower_bound(vref::StructuredVariableRef)
    @assert !JuMP.is_fixed(vref)
    variable_info(vref).lower_bound
end
function JuMP.set_lower_bound(vref::StructuredVariableRef, lower)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(true, lower,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
function JuMP.delete_lower_bound(vref::StructuredVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(false, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.has_upper_bound(vref::StructuredVariableRef) = variable_info(vref).has_ub
function JuMP.upper_bound(vref::StructuredVariableRef)
    @assert !JuMP.is_fixed(vref)
    variable_info(vref).upper_bound
end
function JuMP.set_upper_bound(vref::StructuredVariableRef, upper)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           true, upper,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
function JuMP.delete_upper_bound(vref::StructuredVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           false, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.is_fixed(vref::StructuredVariableRef) = variable_info(vref).has_fix
JuMP.fix_value(vref::StructuredVariableRef) = variable_info(vref).fixed_value
function JuMP.fix(vref::StructuredVariableRef, value)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           true, value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
function JuMP.unfix(vref::StructuredVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           false, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end
JuMP.start_value(vref::StructuredVariableRef) = variable_info(vref).start
function JuMP.set_start_value(vref::StructuredVariableRef, start)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           true, start,
                                           info.binary, info.integer))
end
JuMP.is_binary(vref::StructuredVariableRef) = variable_info(vref).binary
function JuMP.set_binary(vref::StructuredVariableRef)
    @assert !JuMP.is_integer(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           true, info.integer))
end
function JuMP.unset_binary(vref::StructuredVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           false, info.integer))
end
JuMP.is_integer(vref::StructuredVariableRef) = variable_info(vref).integer
function JuMP.set_integer(vref::StructuredVariableRef)
    @assert !JuMP.is_binary(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, true))
end
function JuMP.unset_integer(vref::StructuredVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, false))
end

# Constraints
struct StructuredConstraintRef
    model::StructuredModel # `model` owning the constraint
    idx::Int       # Index in `model.constraints`
end
JuMP.constraint_type(::StructuredModel) = StructuredConstraintRef
function JuMP.add_constraint(m::StructuredModel, c::JuMP.AbstractConstraint, name::String="")
    m.nextconidx += 1
    cref = StructuredConstraintRef(m, m.nextconidx)
    m.constraints[cref.idx] = c
    JuMP.set_name(cref, name)
    cref
end
function MOI.delete!(m::StructuredModel, cref::StructuredConstraintRef)
    delete!(m.constraints, cref.idx)
    delete!(m.connames, cref.idx)
end
MOI.is_valid(m::StructuredModel, cref::StructuredConstraintRef) = cref.idx in keys(m.constraints)
function JuMP.constraint_object(cref::StructuredConstraintRef, F::Type, S::Type)
    c = cref.model.constraints[cref.idx]
    # `TypeError` should be thrown is `F` and `S` are not correct
    # This is needed for the tests in `constraints.jl`
    c.func::F
    c.set::S
    c
end

# Objective

function JuMP.set_objective_sense(m::StructuredModel, sense::MOI.OptimizationSense)
    m.objective_sense = sense
end

function JuMP.set_objective_function(m::StructuredModel, f::AbstractJuMPScalar)
    m.objective_function = f
end

JuMP.set_objective_function(m::StructuredModel, f::Real) = JuMP.set_objective_function(m, convert(AffExpr, f))

JuMP.objective_sense(m::StructuredModel) = m.objective_sense
JuMP.objective_function_type(model::StructuredModel) = typeof(model.objective_function)
JuMP.objective_function(model::StructuredModel) = model.objective_function
function JuMP.objective_function(model::StructuredModel, FT::Type)
    model.objective_function isa FT || error("The objective function is not of type $FT")
    model.objective_function
end

# Names
JuMP.name(vref::StructuredVariableRef) = vref.model.varnames[vref.idx]
function JuMP.set_name(vref::StructuredVariableRef, name::String)
    vref.model.varnames[vref.idx] = name
end
JuMP.name(cref::StructuredConstraintRef) = cref.model.connames[cref.idx]
function JuMP.set_name(cref::StructuredConstraintRef, name::String)
    cref.model.connames[cref.idx] = name
end

# Show
function JuMP.show_backend_summary(io::IO, model::StructuredModel) end
function JuMP.show_objective_function_summary(io::IO, model::StructuredModel)
    println(io, "Objective function type: ",
            JuMP.objective_function_type(model))
end
function JuMP.objective_function_string(print_mode, model::StructuredModel)
    return JuMP.function_string(print_mode, JuMP.objective_function(model))
end
_plural(n) = (isone(n) ? "" : "s")
function JuMP.show_constraints_summary(io::IO, model::StructuredModel)
    n = length(model.constraints)
    print(io, "Constraint", _plural(n), ": ", n)
end
function JuMP.constraints_string(print_mode, model::StructuredModel)
    strings = String[]
    # Sort by creation order, i.e. ConstraintIndex value
    constraints = sort(collect(model.constraints), by = c -> c.first.value)
    for (index, constraint) in constraints
        push!(strings, JuMP.constraint_string(print_mode, constraint))
    end
    return strings
end

include("BendersBridge.jl")

end
