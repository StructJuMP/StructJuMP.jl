import MPI
using JuMP, StochJuMP, DataFrames, Distributions

MPI.init()

N = 100
numScen = 1

m = StochasticModel()

@defVar(m, x[1:N] >= 0)
@addConstraint(m, sum{x[i], i=1:N} >= -1)
@setObjective(m, Min, sum{x[i], i=1:N})

bl = StochasticBlock(m)
@defVar(m, y[1:N] >= 0)
@addConstraint(bl, r=1:N, x[r] + x[r] >= 1/2)
@setObjective(bl, Min, sum{y[i], i=1:N})

StochJuMP.pips_solve(m)

# print(m)
# solve(m)
