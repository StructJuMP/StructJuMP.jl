using StructJuMP

m = StructuredModel()

n = 10
@defVar(m, x[1:n] >= 0)
@defVar(m, y[1:n] >= 0)
for i in 1:n
    @addConstraint(m, x[i] - y[i] <= i)
end
setObjective(m, :Min, x[1] - y[3] + 2x[6] + x[1]*x[2] - 2y[4]*x[3])

numScen = 8

for i in 1:numScen
    bl = StructuredModel(parent=m)
    @defVar(bl, w[1:n] >= 0)
    @addConstraint(bl, sum{i*w[i], i=1:n; iseven(i)} == 1)
    @addConstraint(bl, sum{x[i] - 2w[n-i+1], i=1:n} <= 0)
    setObjective(bl, :Min, w[5] - w[1]*w[1] + 2w[2]*w[3])
end
