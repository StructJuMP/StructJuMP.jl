using StructJuMP, JuMP
using StructJuMPSolverInterface

include("select_solver.jl")

#############
# A sample model
#############
scen = 1
m = StructuredModel(num_scenarios=scen)
@variable(m, x[1:4])
# @variable(m, -100<=x[1:4]<=100)
@NLobjective(m, Min,  1*x[1] + 3*x[3] + 4*x[4] )

for i in 1:scen
    bl = StructuredModel(parent=m)
    # @variable(bl, -200<=y1[1:3]<=200)
    @variable(bl, y1[1:3])
    
    @NLconstraint(bl, -100<= 0.1*x[1] + 0.2*x[2]                          <= 100)
    @NLconstraint(bl, -101<=               1.1*x[2] + 1.2*x[3]            <= 101)
    @NLconstraint(bl, -102<=                          2.1*x[3] + 2.2*x[4] <= 102)

    # @NLconstraint(bl,  9.1*x[1]           + 9.2*x[3]            -109 == 0)
    # @NLconstraint(bl,            8.1*x[2] + 8.2*x[3]            -108 == 0)
    # @NLconstraint(bl,                       7.1*x[3] + 7.2*x[4] -107 == 0)    
    @NLconstraint(bl, 1*x[3] + 2*x[4]  == 107 )

    @NLconstraint(bl, -201<= 10.1*x[1] + 10.2*x[2]                           - 10.3*y1[1] + 10.4*y1[2]               <= 201)
    @NLconstraint(bl, -202<=               11.1*x[2] + 11.2*x[3]             - 11.3*y1[1] + 11.4*y1[2] + 11.5*y1[3]  <= 202)
    @NLconstraint(bl, -203<=                           12.1*x[3] + 12.2*x[4]              + 12.4*y1[2] - 12.5*y1[3]  <= 203)

    @NLobjective(bl, Min, 12*y1[2] + 13*y1[3])
end

structJuMPSolve(m)

getVarValue(m)