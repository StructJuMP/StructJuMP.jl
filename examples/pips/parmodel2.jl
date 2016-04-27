include("../../src/pips_structure_interface.jl")
#an example model

using ParPipsInterface
using StructJuMP, JuMP

scen = 2
firststage = StructuredModel(num_scenarios=scen)
@defVar(firststage, x[1:2])
@addNLConstraint(firststage, x[1] + x[2] == 100)
@setNLObjective(firststage, Min, x[1]^2 + x[2]^2)

for i in 1:scen
    bl = StructuredModel(parent=firststage)
    @defVar(bl, y[1:2])
    @addNLConstraint(bl, x[1] + y[1]+y[2] ≥  0)
    @addNLConstraint(bl, x[1] + y[1]+y[2]  ≤ 50)
    @setNLObjective(bl, Min, y[1]^2 + y[2]^2)
end

ParPipsInterface.solve(firststage)
