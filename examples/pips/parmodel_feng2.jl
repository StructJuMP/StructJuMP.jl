include("../../src/pips_structure_interface.jl")
#an example model

using ParPipsInterface

using StochJuMP, JuMP

#############
# A sample model
#############
scen = 2
m = StochasticModel(num_scenarios=scen)
@defVar(m, x[1:2])
@defVar(m, z[1:2])
@setNLObjective(m, Min, x[1]^2+x[2]^2+x[1]*z[2])
@addNLConstraint(m, x[1] * x[2] <= 100)

for i in 1:scen
    bl = StochasticModel(parent=m)
    @defVar(bl, y[1:2])
    @addNLConstraint(bl, x[2]*x[1] + x[1]*y[1]*y[2] == 10)
    @addNLConstraint(bl, x[2]*x[1]*y[1]*y[2] <= 100)
    @setNLObjective(bl, Min, (y[1]+y[2])^2)
end

ParPipsInterface.solve(m)
