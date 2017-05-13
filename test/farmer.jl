using Base.Test

using Clp

@testset "farmer" begin
    include("../examples/farmer.jl")
    status, objval, soln = DLP(m, ClpSolver())
    @test status == :Optimal
    @test objval == -108390
    @test soln == [170, 80, 250]
end
