# include("../../src/pips_structure_interface.jl")
include("../../src/ipopt_interface.jl")
#an example model

# using ParPipsInterface
using SerialIpoptInterface

using StructJuMP, JuMP

#############
# A sample model
#############
scen = 2
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:2])
@setNLObjective(m, Min, x[1]^2+x[2]^2)
@addNLConstraint(m, x[1] * x[2] <= 100)

for i in 1:scen
    bl = StructuredModel(parent=m)
    @defVar(bl, y)
    @addNLConstraint(bl, x[2]*x[1] + x[1]*y <= 10)
    @setNLObjective(bl, Min, y^2)
end

# ParPipsInterface.structJuMPSolve(m)
SerialIpoptInterface.structJuMPSolve(m)

getVarValue(m)

