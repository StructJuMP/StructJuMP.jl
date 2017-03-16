using StructJuMP

firststage = StructuredModel()
@variable(firststage, x[1:2])
@constraint(firststage, sum(x) == 100)
@NLobjective(firststage, :Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for scen in 1:100
    bl = StructuredModel(parent=firststage, id=scen)
    @variable(bl, y[1:2])
    idx = (isodd(scen) ? 1 : 2)
    @constraint(bl, x[idx] + sum(y) ≥  0)
    @constraint(bl, x[idx] + sum(y) ≤ 50)
    @NLobjective(bl, :Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end
