using JuMP
using StructJuMP
using ECOS
using FactCheck
using Cbc

include("../src/BendersBridge.jl")

misocp_solver = CbcSolver()
socp_solver = ECOS.ECOSSolver()

facts("[Benders] Empty scenario test") do

    m = StructuredModel(num_scenarios=0)
    @variable(m, x, Int)
    @constraint(m, x <= 4)
    @objective(m, :Min, -5*x)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-20.0)

end

facts("[Benders] Infeasible problem test") do

    numScen = 1
    m = StructuredModel(num_scenarios=numScen)

    @variable(m, x, Int)

    @constraint(m, x <= 1)
    @objective(m, :Min, -5*x)

    bl = StructuredModel(parent=m, id=1)
    @variable(bl, y1 >= 2)
    @variable(bl, y2 <= 2)
    @constraint(bl, x >= y1)
    @constraint(bl, norm(y1) <= y2)
    @objective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Infeasible

end

facts("[Benders] Infeasibility cut execution test #1") do

    numScen = 1
    m = StructuredModel(num_scenarios=numScen)

    @variable(m, x, Int)

    @constraint(m, x <= 4)
    @objective(m, :Min, -5*x)

    bl = StructuredModel(parent=m, id=1)
    @variable(bl, y1 >= 0)
    @variable(bl, y2 <= 2)
    @constraint(bl, x <= y1)
    @constraint(bl, norm(y1) <= y2)
    @objective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-4.0)

end

facts("[Benders] Optimality cut execution test #1") do

    numScen = 1
    m = StructuredModel(num_scenarios=numScen)

    @variable(m, x, Int)

    @constraint(m, x <= 4)
    @objective(m, :Min, -5*x)

    bl = StructuredModel(parent=m, id=1)
    @variable(bl, y1 >= 2)
    @variable(bl, y2 <= 4)
    @constraint(bl, x <= y1)
    @constraint(bl, norm(y1) <= y2)
    @objective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-8.0)

end
