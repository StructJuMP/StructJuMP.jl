using Test

using StructJuMP

@testset "printhook" begin

    numScen = 2
    parent = StructuredModel(num_scenarios=numScen)

    @variable(parent, 0 <= x <= 1)
    @variable(parent, 0 <= y <= 1)

    @constraint(parent, x + y == 1)
    @objective(parent, Min, x*x + y)

    for i in 1:numScen
        bl = StructuredModel(parent=parent, id=i)
        @variable(bl, w >= 0)
        @constraint(bl, w - x - y <= 1)
        @objective(bl, Min, w*w + w)
    end

    str = string(parent)
    @test sprint(print, parent) == """
Min x² + y
Subject to
 x + y = 1.0
"""
    @test sprint(show, parent) == """
A JuMP Model
Minimization problem with:
Variables: 2
Objective function type: GenericQuadExpr{Float64,StructJuMP.StructuredVariableRef}
Constraint: 1
Names registered in the model: x, y"""
    @test sprint(show, "text/latex", parent) == """
\$\$ \\begin{alignat*}{1}\\min\\quad & x^2 + y\\\\
\\text{Subject to} \\quad & x + y = 1.0\\\\
\\end{alignat*}
 \$\$"""
    #@test occursin("Child", str)
end
