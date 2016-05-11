include("../../src/pips_structure_interface.jl")
# include("../../src/ipopt_interface.jl")
#an example model

using ParPipsInterface
# using SerialIpoptInterface

using StructJuMP, JuMP

#############
# A sample model
#############
scen = 1
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:4])
# @defVar(m, -100<=x[1:4]<=100)
@setNLObjective(m, Min,  1*x[1] + 3*x[3] + 4*x[4] )

for i in 1:scen
    bl = StructuredModel(parent=m)
    # @defVar(bl, -200<=y1[1:3]<=200)
    @defVar(bl, y1[1:3])
    
    @addNLConstraint(bl, -100<= 0.1*x[1] + 0.2*x[2]                          <= 100)
    @addNLConstraint(bl, -101<=               1.1*x[2] + 1.2*x[3]            <= 101)
    @addNLConstraint(bl, -102<=                          2.1*x[3] + 2.2*x[4] <= 102)

    # @addNLConstraint(bl,  9.1*x[1]           + 9.2*x[3]            -109 == 0)
    # @addNLConstraint(bl,            8.1*x[2] + 8.2*x[3]            -108 == 0)
    # @addNLConstraint(bl,                       7.1*x[3] + 7.2*x[4] -107 == 0)    
    @addNLConstraint(bl, 1*x[3] + 2*x[4]  == 107 )

    @addNLConstraint(bl, -201<= 10.1*x[1] + 10.2*x[2]                           - 10.3*y1[1] + 10.4*y1[2]               <= 201)
    @addNLConstraint(bl, -202<=               11.1*x[2] + 11.2*x[3]             - 11.3*y1[1] + 11.4*y1[2] + 11.5*y1[3]  <= 202)
    @addNLConstraint(bl, -203<=                           12.1*x[3] + 12.2*x[4]              + 12.4*y1[2] - 12.5*y1[3]  <= 203)

    @setNLObjective(bl, Min, 12*y1[2] + 13*y1[3])
end

ParPipsInterface.structJuMPSolve(m)
# SerialIpoptInterface.structJuMPSolve(m)

getVarValue(m)