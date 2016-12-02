using StructJuMP

N = 1000
numScen = 1

m = StructuredModel(numScen)

@variable(m, x[1:N] >= 0)
@constraint(m, sum{x[i], i=1:N} >= -1)
@JuMP.setobjective(m, Min, sum{x[i], i=1:N})

bl = StructuredModel(parent=m)
@variable(bl, y[1:N] >= 0)
@constraint(bl, constr[r=2:(N-1)], x[r] + x[r-1] + y[r] + y[r+1] >= 0)
@JuMP.setobjective(bl, Min, sum{y[i], i=1:N})
