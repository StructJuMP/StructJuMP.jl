using StructJuMP, JuMP
using StructJuMPSolverInterface

include("select_solver.jl")
#an example model
scen = 2
m = StructuredModel(num_scenarios=scen)
@variable(m, x[1:2])
@NLobjective(m, Min, (x[1]+x[2])^2)
@NLconstraint(m, x[1] * x[2] == 10)

for i in 1:scen
    bl = StructuredModel(parent=m)
    @variable(bl, y)
    @NLconstraint(bl, x[2]^2 + x[1]*y ≤ 5)
    @NLobjective(bl, Min, (x[1]+x[2])*y)
end

structJuMPSolve(m)

getVarValue(m)