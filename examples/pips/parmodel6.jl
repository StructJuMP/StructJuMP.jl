using StructJuMP, JuMP
using SolverInterface

include("select_solver.jl")

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
    @defVar(bl, y[1:2])
    @addNLConstraint(bl, x[2]*x[1] + x[1]*y[1]*y[2] <= 10)
    @setNLObjective(bl, Min, (y[1]+y[2])^2)
end

structJuMPSolve(m)

getVarValue(m)
