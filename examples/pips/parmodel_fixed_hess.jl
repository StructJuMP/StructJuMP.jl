using StructJuMP, JuMP
using StructJuMPSolverInterface

#############
# A sample model
#############
scen = 1
m = StructuredModel(num_scenarios=scen)
@variable(m, -100<=x[1:4]<=100)
@NLobjective(m, Min,  0.1*x[1]^2 + 0.2*x[2]^2 + 0.3*x[3]^2 + 0.4*x[4]^2 +
    + 1.1*x[1]*x[2] + 2.1*x[1]*x[3] + 2.2*x[1]*x[4] + 
    + 2.1*x[2]*x[3] + 2.2*x[2]*x[4] 
    + 3.1*x[3]*x[4])

for i in getLocalChildrenIds(m)
    bl = StructuredModel(parent=m, id=i)
    @variable(bl, -200<=y1[1:3]<=200)
    @NLconstraint(bl, -100<= 0.1*x[1] + 0.2*x[2]                          <= 100)
    @NLconstraint(bl, -101<=               1.1*x[2] + 1.2*x[3]            <= 101)
    @NLconstraint(bl, -102<=                          2.1*x[3] + 2.2*x[4] <= 102)


    @NLconstraint(bl, -201<= 10.1*x[1] + 10.2*x[2]                           - 10.3*y1[1] + 10.4*y1[2]               <= 201)
    @NLconstraint(bl, -202<=               11.1*x[2] + 11.2*x[3]             - 11.3*y1[1] + 11.4*y1[2] + 11.5*y1[3]  <= 202)
    @NLconstraint(bl, -203<=                           12.1*x[3] + 12.2*x[4]              + 12.4*y1[2] - 12.5*y1[3]  <= 203)

    @NLobjective(bl, Min, 
        5.1*y1[1]*x[1] + 5.2*y1[1]*x[2] + 5.3*y1[1]*x[3] + 5.4*y1[1]*x[4]
    +                    6.2*y1[2]*x[2] +                  6.4*y1[2]*x[4]
    + 7.1*y1[3]*x[1] +                  7.3*y1[3]*x[3]
    + 8.1*y1[1]^2    + 8.2*y1[2]^2 + 8.3*y1[3]^2
    + 9.1*y1[1]*y1[2]+ 9.2*y1[1]*y1[3]
    + 10.1*y1[2]*y1[3])
end

solve(m, solver="PipsNlp")
