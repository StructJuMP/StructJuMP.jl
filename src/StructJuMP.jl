module StructJuMP
using JuMP # To reexport, should be using (not import)
using MathOptInterface
const MOI = MathOptInterface

mutable struct StructuredModel <: AbstractModel
    # special variablewise properties that we keep track of:
    # lower bound, upper bound, fixed, integrality, binary
    variabletolowerbound::Dict{MOIVAR,MOILB}
    variabletoupperbound::Dict{MOIVAR,MOIUB}
    variabletofix::Dict{MOIVAR,MOIFIX}
    variabletointegrality::Dict{MOIVAR,MOIINT}
    variabletozeroone::Dict{MOIVAR,MOIBIN}

    customnames::Vector

    objdict::Dict{Symbol,Any} # dictionary from variable and constraint names to objects
end

end
