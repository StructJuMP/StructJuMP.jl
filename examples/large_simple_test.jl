using StructJuMP

N = 1000
numScen = 1

m = StructuredModel(num_scenarios=numScen)

@variable(m, x[1:N] >= 0)
@constraint(m, sum(x[i] for i=1:N) >= -1)
@objective(m, :Min, sum(x[i] for i=1:N))

bl = StructuredModel(parent=m, id=1)
@variable(bl, y[1:N] >= 0)
@constraint(bl, constr[r=2:(N-1)], x[r] + x[r-1] + y[r] + y[r+1] >= 0)
@objective(bl, :Min, sum(y[i] for i=1:N))
