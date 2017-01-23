using StructJuMP
using Base.Test


@testset "printhook" begin

    numScen = 2
    m = StructuredModel(num_scenarios=numScen)

    @variable(m, 0 <= x <= 1)
    @variable(m, 0 <= y <= 1)

    @constraint(m, x + y == 1)
    @objective(m, :Min, x*x + y)

    for i in 1:numScen
        bl = StructuredModel(parent=m, id=i)
        @variable(bl, w >= 0)
        @constraint(bl, w - x - y <= 1)
        @objective(bl, :Min, w*w + w)
    end

    str = string(m)
    @test contains(str, "Child")
end
