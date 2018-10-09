using Compat
using Compat.Test

using GLPKMathProgInterface

@testset "farmer" begin
    include("../examples/farmer.jl")
    status, objval, soln = DLP(m, GLPKSolverLP())
    @test status == :Optimal
    @test objval ≈ -108390
    @test soln ≈ [170, 80, 250]
end
