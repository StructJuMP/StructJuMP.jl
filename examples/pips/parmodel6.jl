using StructJuMP, JuMP
using StructJuMPSolverInterface

include("select_solver.jl")

#############
# A sample model
#############
scen = 2
m = StructuredModel(num_scenarios=scen)
@variable(m, x[1:2])
@NLobjective(m, Min, x[1]^2+x[2]^2)
@NLconstraint(m, x[1] * x[2] <= 100)

for i in 1:scen
    bl = StructuredModel(parent=m)
    @variable(bl, y[1:2])
    @NLconstraint(bl, x[2]*x[1] + x[1]*y[1]*y[2] <= 10)
    @NLobjective(bl, Min, (y[1]+y[2])^2)
end

structJuMPSolve(m)

getVarValue(m)
