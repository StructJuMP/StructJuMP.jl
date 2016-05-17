using StructJuMP, JuMP
using SolverInterface

include("select_solver.jl")

#an example model

nx1=50; ni1=10
nx2=5000; ni2=30

scen = 30
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:nx1])
@setNLObjective(m, Min, sum{(x[i]-1)^4, i=1:nx1})

#@addNLConstraint(m, x[1] + x[2] >= 3)


for i in 1:scen
    bl = StructuredModel(parent=m)
    @defVar(bl, y[1:nx2])
    #@addNLConstraint(bl, x[1] + y[1]+y[2] >= -10)
    #@addNLConstraint(bl, x[2] + y[1]+y[2] >= -50)
    @setNLObjective(bl, Min, sum{y[k]^4, k=1:nx2})

    for c in 1:ni2
      @addNLConstraint(bl,(x[(c-1)%nx1 + 1] + y[(c-1)%nx2+1])^2 <= 16)
    end
end

structJuMPSolve(m)

getVarValue(m)