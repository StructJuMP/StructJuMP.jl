using Test


include("../examples/easy_test.jl")

@testset "easy_test" begin
    include("../examples/easy_test.jl")
end

@testset "harder_test" begin
    include("../examples/harder_test.jl")
end

@testset "large_simple_test" begin
    include("../examples/large_simple_test.jl")
end

@testset "explicit_stochastic" begin
    include("../examples/explicit_stochastic.jl")
end

@testset "transportation" begin
    include("../examples/transportation.jl")
end

@testset "parmodel1" begin
    include("../examples/parmodel1.jl")
    # using Ipopt
    # set_optimizer(firststage, Ipopt.Optimizer)
    # set_optimizer_attribute(first_stage, "print_level", 0)
    # optimize!(firststage)
end

@testset "parmodel2" begin
    include("../examples/parmodel2.jl")
    # using Ipopt
    # set_optimizer(firststage, Ipopt.Optimizer)
    # set_optimizer_attribute(first_stage, "print_level", 0)
    # optimize!(firststage)
end

@testset "parmodel3" begin
    include("../examples/parmodel3.jl")
    # using Ipopt
    # set_optimizer(firststage, Ipopt.Optimizer)
    # set_optimizer_attribute(first_stage, "print_level", 0)
    # optimize!(firststage)
end

@testset "parmodel4" begin
    include("../examples/parmodel4.jl")
    # using Ipopt
    # set_optimizer(firststage, Ipopt.Optimizer)
    # set_optimizer_attribute(first_stage, "print_level", 0)
    # optimize!(firststage)
end

@testset "farmer" begin
    include("../examples/farmer.jl")
    # status, objval, soln = DLP(m, GLPK.Optimizer)
    # @test status == :Optimal
    # @test objval ≈ -108390
    # @test soln ≈ [170, 80, 250]
end

#=
@testset "stochpdegas" begin
    include("../examples/stochpdegas.jl")
end
=#

#=
@testset "simple" begin
    include("../examples/simple.jl")
end
=#
