using JuMP
using StochJuMP
using CPLEX
using ECOS
using FactCheck

include("../src/BendersBridge.jl")

misocp_solver = CplexSolver()
socp_solver = ECOS.ECOSSolver()

facts("[Benders] Empty scenario test") do

    m = StochasticModel(0)
    @defVar(m, x, Int)
    @addConstraint(m, x <= 4)
    @setObjective(m, :Min, -5*x)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-20.0)

end

facts("[Benders] Infeasible problem test") do

    numScen = 1
    m = StochasticModel(numScen)

    @defVar(m, x, Int)

    @addConstraint(m, x <= 1)
    @setObjective(m, :Min, -5*x)

    bl = StochasticModel(parent=m)
    @defVar(bl, y1 >= 2)
    @defVar(bl, y2 <= 2)
    @addConstraint(bl, x >= y1)
    @addConstraint(bl, norm(y1) <= y2)
    @setObjective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Infeasible

end

facts("[Benders] Infeasibility cut execution test #1") do

    numScen = 1
    m = StochasticModel(numScen)

    @defVar(m, x, Int)

    @addConstraint(m, x <= 4)
    @setObjective(m, :Min, -5*x)

    bl = StochasticModel(parent=m)
    @defVar(bl, y1 >= 0)
    @defVar(bl, y2 <= 2)
    @addConstraint(bl, x <= y1)
    @addConstraint(bl, norm(y1) <= y2)
    @setObjective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-4.0)

end

facts("[Benders] Optimality cut execution test #1") do

    numScen = 1
    m = StochasticModel(numScen)

    @defVar(m, x, Int)

    @addConstraint(m, x <= 4)
    @setObjective(m, :Min, -5*x)

    bl = StochasticModel(parent=m)
    @defVar(bl, y1 >= 2)
    @defVar(bl, y2 <= 4)
    @addConstraint(bl, x <= y1)
    @addConstraint(bl, norm(y1) <= y2)
    @setObjective(bl, Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @fact output[1] --> :Optimal
    @fact output[2] --> roughly(-8.0)

end
