using StructJuMP

N = 1000
numScen = 1

m = StructuredModel(numScen)

@defVar(m, x[1:N] >= 0)
@addConstraint(m, sum{x[i], i=1:N} >= -1)
@setObjective(m, Min, sum{x[i], i=1:N})

bl = StructuredModel(parent=m)
@defVar(bl, y[1:N] >= 0)
@addConstraint(bl, constr[r=2:(N-1)], x[r] + x[r-1] + y[r] + y[r+1] >= 0)
@setObjective(bl, Min, sum{y[i], i=1:N})
