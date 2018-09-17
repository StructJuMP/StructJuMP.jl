module StructJuMP

export StructuredModel

using Compat

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
    parent::Union{StructuredModel, Compat.Nothing}
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
    objective_sense::Symbol
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
                    :Min, zero(JuMP.GenericAffExpr{Float64, StructuredVariableRef}), # Model objective
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
default_probability(::Compat.Nothing) = 1.0

#### Model ####

JuMP.object_dictionary(m::StructuredModel) = m.objdict

# Variables
struct StructuredVariableRef <: JuMP.AbstractVariableRef
    model::StructuredModel # `model` owning the variable
    idx::Int       # Index in `model.variables`
end
Base.copy(v::StructuredVariableRef) = v
Base.:(==)(v::StructuredVariableRef, w::StructuredVariableRef) = v.model === w.model && v.idx == w.idx
JuMP.owner_model(v::StructuredVariableRef) = v.model
JuMP.isequal_canonical(v::StructuredVariableRef, w::StructuredVariableRef) = v == w
JuMP.variabletype(::StructuredModel) = StructuredVariableRef
function JuMP.add_variable(m::StructuredModel, v::JuMP.AbstractVariable, name::String="")
    m.nextvaridx += 1
    vref = StructuredVariableRef(m, m.nextvaridx)
    m.variables[vref.idx] = v
    JuMP.setname(vref, name)
    vref
end
function MOI.delete!(m::StructuredModel, vref::StructuredVariableRef)
    delete!(m.variables, vref.idx)
    delete!(m.varnames, vref.idx)
end
MOI.is_valid(m::StructuredModel, vref::StructuredVariableRef) = vref.idx in keys(m.variables)
JuMP.num_variables(m::StructuredModel) = length(m.variables)

JuMP.has_lower_bound(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.haslb
function JuMP.lower_bound(vref::StructuredVariableRef)
    @assert !JuMP.is_fixed(vref)
    vref.model.variables[vref.idx].info.lowerbound
end
function JuMP.set_lower_bound(vref::StructuredVariableRef, lower)
    vref.model.variables[vref.idx].info.haslb = true
    vref.model.variables[vref.idx].info.lowerbound = lower
end
JuMP.has_upper_bound(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.hasub
function JuMP.upper_bound(vref::StructuredVariableRef)
    @assert !JuMP.isfixed(vref)
    vref.model.variables[vref.idx].info.upperbound
end
function JuMP.set_upper_bound(vref::StructuredVariableRef, upper)
    vref.model.variables[vref.idx].info.hasub = true
    vref.model.variables[vref.idx].info.upperbound = upper
end
JuMP.is_fixed(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.hasfix
JuMP.fix_value(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.fixedvalue
function JuMP.fix(vref::StructuredVariableRef, value)
    vref.model.variables[vref.idx].info.fixedvalue = value
end
JuMP.start_value(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.start
function JuMP.set_start_value(vref::StructuredVariableRef, start)
    vref.model.variables[vref.idx].info.start = start
end
JuMP.is_binary(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.binary
function JuMP.set_binary(vref::StructuredVariableRef)
    @assert !JuMP.isinteger(vref)
    vref.model.variables[vref.idx].info.binary = true
end
function JuMP.unset_binary(vref::StructuredVariableRef)
    vref.model.variables[vref.idx].info.binary = false
end
JuMP.is_integer(vref::StructuredVariableRef) = vref.model.variables[vref.idx].info.integer
function JuMP.set_integer(vref::StructuredVariableRef)
    @assert !JuMP.isbinary(vref)
    vref.model.variables[vref.idx].info.integer = true
end
function JuMP.unset_integer(vref::StructuredVariableRef)
    vref.model.variables[vref.idx].info.integer = false
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
    JuMP.setname(cref, name)
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
function JuMP.set_objective(m::StructuredModel, sense::Symbol, f::JuMP.AbstractJuMPScalar)
    m.objective_sense = sense
    m.objective_function = f
end
JuMP.objective_sense(m::StructuredModel) = m.objective_sense
function JuMP.objective_function(m::StructuredModel, FT::Type)
    m.objective_function isa FT || error("The objective function is not of type $FT")
    m.objective_function
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

include("BendersBridge.jl")

end
