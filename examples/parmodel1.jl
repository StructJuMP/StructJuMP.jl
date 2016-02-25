using StochJuMP, JuMP

firststage = StochasticModel()
@defVar(firststage, x[1:2])
@setNLObjective(firststage, Min, (x[1]+x[2])^2)
@addNLConstraint(firststage, x[1] * x[2] == 10)

for scen in 1:2
    bl = StochasticModel(parent=firststage)
    @defVar(bl, y)
    @addConstraint(bl, x[2]^2 + x[1]*y ≤ 5)
    @setNLObjective(bl, Min, (x[1]+x[2])*y)
end

solve(firststage)
