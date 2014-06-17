using JuMP

m = Model()

@defVar(m, x >= 1)
@defVar(m, y <= 2)

@addConstraint(m,  x + y == 1)
@addConstraint(m, -x + y <= 1)

@defVar(m, 0 <= w[1:2] <= 1)
@addConstraint(m, w[1] - x - y <= 1)
@addConstraint(m, w[1] + x == 2)

# @addConstraint(m, 2w[2] - x - y <= 5)
# @addConstraint(m, 3w[2] + x == 4)

# setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y + w[1]*w[1] + w[1] + 0.5w[2]*w[2] + 0.75w[2])
setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y + w[1]*w[1] + w[1])

solve(m)
