using StructJuMP

firststage = StructuredModel()
@variable(firststage, x[1:2])
@objective(firststage, Min, (x[1]+x[2])^2)
@constraint(firststage, x[1] * x[2] == 10)

for scen in 1:2
    local bl = StructuredModel(parent=firststage, id=scen)
    @variable(bl, y)
    @constraint(bl, x[2]^2 + x[1]*y ≤ 5)
    @objective(bl, Min, (x[1]+x[2])*y)
end
