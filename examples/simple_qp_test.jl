using StructJuMP

m = StructuredModel()

n = 10
@variable(m, x[1:n] >= 0)
@variable(m, y[1:n] >= 0)
for i in 1:n
    @constraint(m, x[i] - y[i] <= i)
end
@objective(m, :Min, x[1] - y[3] + 2x[6] + x[1]*x[2] - 2y[4]*x[3])

numScen = 8

for i in 1:numScen
    bl = StructuredModel(parent=m,id=i)
    @variable(bl, w[1:n] >= 0)
    @constraint(bl, sum(i*w[i] for i=2:2:n) == 1)
    @constraint(bl, sum(x[i] - 2w[n-i+1] for i=1:n) <= 0)
    @objective(bl, :Min, w[5] - w[1]*w[1] + 2w[2]*w[3])
end
