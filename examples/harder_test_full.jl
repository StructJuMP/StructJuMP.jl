using JuMP

m = Model()

@defVar(m, x >= 1)
@defVar(m, y <= 2)

@addConstraint(m,  x + y == 1)
@addConstraint(m, -x + y <= 1)

rhs  = [5,4]
coef = [2,3]
qc   = [1,0.5]
ac   = [1,0.75]

@defVar(m, 0 <= w[1:2] <= 1)
@addConstraint(m, coef[1]*w[1] - x - y <= rhs[1])
@addConstraint(m, coef[1]*w[1] + x     == rhs[1])

@addConstraint(m, coef[2]*w[2] - x - y <= rhs[2])
@addConstraint(m, coef[2]*w[2] + x     == rhs[2])

setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y + 
                      qc[1]*w[1]*w[1] + ac[1]*w[1] + 
                      qc[2]*w[2]*w[2] + ac[2]*0.75w[2])

solve(m)
