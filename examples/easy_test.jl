using StructJuMP

numScen = 2
m = StructuredModel(numScen)

@variable(m, 0 <= x <= 1)
@variable(m, 0 <= y <= 1)

@constraint(m, x + y == 1)
@objective(m, :Min, x*x + y)

for i in 1:numScen
    bl = StructuredModel(parent=m)
    @variable(bl, w >= 0)
    @constraint(bl, w - x - y <= 1)
    @objective(bl, :Min, w*w + w)
end
