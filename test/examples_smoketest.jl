using Compat, Compat.Test


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
    using Ipopt
    setsolver(firststage, IpoptSolver(print_level=0))
    solve(firststage)
end

@testset "parmodel2" begin
    include("../examples/parmodel2.jl")
    using Ipopt
    setsolver(firststage, IpoptSolver(print_level=0))
    solve(firststage)
end

@testset "parmodel3" begin
    include("../examples/parmodel3.jl")
    using Ipopt
    setsolver(firststage, IpoptSolver(print_level=0))
    solve(firststage)
end

@testset "parmodel4" begin
    include("../examples/parmodel4.jl")
    using Ipopt
    setsolver(firststage, IpoptSolver(print_level=0))
    solve(firststage)
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
