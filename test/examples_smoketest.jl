using Base.Test


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

#=
@testset "simple" begin
    include("../examples/simple.jl")
end
=#
