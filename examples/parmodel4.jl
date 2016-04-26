include("../src/pips_structure_interface.jl")
#an example model

using ParPipsInterface
using StructJuMP, JuMP

scen = 10
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:2])
@addConstraint(m, sum(x) == 100)
@setNLObjective(m, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for i in 1:scen
    bl = StructuredModel(parent=m)
    @defVar(bl, y[1:2])
    idx = (isodd(i) ? 1 : 2)
    @addNLConstraint(bl, x[idx] + y[1]+y[2] ≥  0)
    @addNLConstraint(bl, x[idx] + y[1]+y[2] ≤ 50)
    @setNLObjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end


ParPipsInterface.solve(m)
@show "done"

