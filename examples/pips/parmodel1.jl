include("../../src/pips_structure_interface.jl")
#an example model

using ParPipsInterface
using StructJuMP, JuMP

scen = 2
firststage = StructuredModel(num_scenarios=scen)
@defVar(firststage, x[1:2])
@setNLObjective(firststage, Min, (x[1]+x[2])^2)
@addNLConstraint(firststage, x[1] * x[2] == 10)

for i in 1:scen
    bl = StructuredModel(parent=firststage)
    @defVar(bl, y)
    @addNLConstraint(bl, x[2]^2 + x[1]*y ≤ 5)
    @setNLObjective(bl, Min, (x[1]+x[2])*y)
end

ParPipsInterface.solve(firststage)
