import MPI # put this first!
using StochJuMP, JuMP

MPI.init()

m = StochasticModel()

@defVar(m, x >= 1)
@defVar(m, y <= 2)

@addConstraint(m,  x + y == 1)
@addConstraint(m, -x + y <= 1)
setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y)

numScen = 2

rhs  = [5,4]
coef = [2,3]

for i in 1:numScen
    bl = StochasticBlock(m, numScen)
    @defVar(bl, 0 <= w <= 1)
    @addConstraint(bl, coef[i]w - x - y <= rhs)
    @addConstraint(bl, coef[i]w + x     == rhs)
    setObjective(bl, :Min, w*w + w)
end

StochJuMP.pips_solve(m)
