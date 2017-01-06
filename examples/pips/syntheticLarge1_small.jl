using StructJuMP, JuMP
using StructJuMPSolverInterface

include("select_solver.jl")

#an example model
nx1=5; ni1=3
nx2=3; ni2=3

scen = 2
m = StructuredModel(num_scenarios=scen)
@variable(m, x[1:nx1])
@NLobjective(m, Min, sum((x[i]-1)^4 for i=1:nx1))

#@NLconstraint(m, x[1] + x[2] >= 3)


for i in 1:scen
    bl = StructuredModel(parent=m)
    @variable(bl, y[1:nx2])
    #@NLconstraint(bl, x[1] + y[1]+y[2] >= -10)
    #@NLconstraint(bl, x[2] + y[1]+y[2] >= -50)
    @NLobjective(bl, Min, sum(y[k]^4 for k=1:nx2))

    for c in 1:ni2
      @NLconstraint(bl,(x[(c-1)%nx1 + 1] + y[(c-1)%nx2+1])^2 <= 16)
    end
end

structJuMPSolve(m)

getVarValue(m)

# #verification 
# nx1=5
# nx2=3
# sum((x-1).^4) + sum(y1.^4) + sum(y2.^4)

# for c in 1:3
# @show (x[(c-1)%5 + 1] + y1[(c-1)%3+1])^2
# end
# for c in 1:3
# @show (x[(c-1)%5 + 1] + y2[(c-1)%3+1])^2
# end