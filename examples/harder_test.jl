using StructJuMP

numScens = 2
m = StructuredModel(num_scenarios=numScens)

@variable(m, x >= 1)
@variable(m, y <= 2)

@constraint(m,  x + y == 1)
@constraint(m, -x + y <= 1)
@objective(m, :Min, x*x + 0.5x*y + 0.25y*y - y)

rhs  = [5,4]
coef = [2,3]
qc   = [1,0.5]
ac   = [1,0.75]

for i = 1:numScens
    bl = StructuredModel(parent=m, id=i)
    @variable(bl, 0 <= w <= 1)
    @constraint(bl, coef[i]w - x - y <= rhs[i])
    @constraint(bl, coef[i]w + x     == rhs[i])
    @objective(bl, :Min, qc[i]*w*w + ac[i]*w)
end
