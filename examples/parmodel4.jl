using StructJuMP, JuMP

firststage = StructuredModel()
@defVar(firststage, x[1:2])
@addConstraint(firststage, sum(x) == 100)
@setNLObjective(firststage, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for scen in 1:100
    bl = StructuredModel(parent=firststage)
    @defVar(bl, y[1:2])
    idx = (isodd(scen) ? 1 : 2)
    @addConstraint(bl, x[idx] + sum(y) ≥  0)
    @addConstraint(bl, x[idx] + sum(y) ≤ 50)
    @setNLObjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end

solve(firststage)
