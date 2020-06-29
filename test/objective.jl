using Test
using StructJuMP

@testset "objective functions" begin
    m = StructuredModel(num_scenarios = 2)
    @variable(m, x >= 0)
    @objective(m, Min, x*x)

    sm = StructuredModel(parent = m, id = 1)
    @variable(sm, y >= 0)
    @objective(sm, Max, y)

    sm = StructuredModel(parent = m, id = 2)
    @variable(sm, y >= 0)
    @objective(sm, Min, 0)

    @test objective_sense(m) == MOI.MIN_SENSE
    @test objective_sense(m.children[1]) == MOI.MAX_SENSE
    @test objective_sense(m.children[2]) == MOI.MIN_SENSE
    @test objective_function_type(m) == GenericQuadExpr{Float64,StructJuMP.StructuredVariableRef}
    @test objective_function_type(m.children[1]) == StructJuMP.StructuredVariableRef
    @test objective_function_type(m.children[2]) == GenericAffExpr{Float64,VariableRef}
end