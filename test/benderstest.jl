using JuMP
using StructJuMP
using ECOS
using Cbc
using Base.Test

misocp_solver = CbcSolver()
socp_solver = ECOS.ECOSSolver(verbose=false)

@testset "[Benders] conicconstraintdata with more variables in parent" begin
    m = StructuredModel(num_scenarios=1)
    @variable(m, x[1:2])
    @objective(m, :Min, sum(x))

    bl = StructuredModel(parent=m, id=1)
    @variable(bl, y)
    @constraint(bl, 4y + 5x[1] + 6x[2] >= 2)
    @objective(bl, :Max, 3y)

    c, A, B, b, var_cones, con_cones, v = StructJuMP.conicconstraintdata(bl)
    @test c == [-3]
    @test A == [5 6]
    @test B == reshape([4], 1, 1)
    @test b == [2]
    @test length(var_cones) == 1
    @test var_cones[1][1] == :Free
    @test collect(var_cones[1][2]) == [1]
    @test con_cones[1][1] == :NonPos
    @test collect(con_cones[1][2]) == [1]
    @test v == [:Cont]
end

@testset "[Benders] Empty scenario test" begin

    m = StructuredModel(num_scenarios=0)
    @variable(m, x, Int)
    @constraint(m, x <= 4)
    @objective(m, :Min, -5*x)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @test output[1] == :Optimal
    @test isapprox(output[2], -20.0)

end

@testset "[Benders] Infeasible problem test" begin

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
    @objective(bl, :Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @test output[1] == :Infeasible

end

@testset "[Benders] Infeasibility cut execution test #1" begin

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
    @objective(bl, :Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @test output[1] == :Optimal
    @test isapprox(output[2], -4.0)

end

@testset "[Benders] Optimality cut execution test #1" begin

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
    @objective(bl, :Min, 2*y1 + y2)

    output = BendersBridge(m, misocp_solver, socp_solver)

    @test output[1] == :Optimal
    @test isapprox(output[2], -8.0)

end
