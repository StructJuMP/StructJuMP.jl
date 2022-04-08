using StructJuMP

firststage = StructuredModel()
@variable(firststage, x[1:2])
@constraint(firststage, sum(x) == 100)
@objective(firststage, Min, x[1]^2 + x[2]^2)

for scen in 1:2
    local bl = StructuredModel(parent=firststage, id=scen)
    @variable(bl, y[1:2])
    @constraint(bl, x[3-scen] + sum(y) ≥  0)
    @constraint(bl, x[3-scen] + sum(y) ≤ 50)
    @objective(bl, Min, y[1]^2 + y[2]^2)
end
