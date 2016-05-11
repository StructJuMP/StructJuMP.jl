# include("../../src/pips_structure_interface.jl")
include("../../src/ipopt_interface.jl")
#an example model

# using ParPipsInterface
using SerialIpoptInterface

using StructJuMP, JuMP


nx1=5; ni1=3
nx2=3; ni2=3

scen = 2
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

# ParPipsInterface.structJuMPSolve(m)
SerialIpoptInterface.structJuMPSolve(m)

getVarValue(m)