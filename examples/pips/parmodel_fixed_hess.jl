include("../../src/pips_structure_interface.jl")
#an example model

using ParPipsInterface

using StructJuMP, JuMP

#############
# A sample model
#############
scen = 2
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:4])
@setNLObjective(m, Min, x[1]^2+ 2*x[2]^2 + 3*x[3]^2 + 4*x[4]^2 + x[1]*x[2] + x[1]*x[3] + x[1]*x[4] + x[2]*x[3] + x[2]*x[4] + x[3]*x[4])
@addNLConstraint(m, x[1] + x[2]<= 101)
@addNLConstraint(m, x[2] + x[3]<= 102)
@addNLConstraint(m, x[3] + x[4]<= 103)

for i in 1:scen
    bl = StructuredModel(parent=m)
    @defVar(bl, y[1:3])
    @addNLConstraint(bl, x[1] + x[2] - y[1] + y[2] <= 201)
    @addNLConstraint(bl, x[2] + x[3] - y[2] + y[3] <= 202)
    @addNLConstraint(bl, x[3] + x[4] - y[3]        <= 203)

    @setNLObjective(bl, Min, y[1]^2 + 2*y[2]^2 + 3*y[3]^2 + y[1]*y[2] + y[1]*y[3] + y[2]*y[3])
end

ParPipsInterface.solve(m)
