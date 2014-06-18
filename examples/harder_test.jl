import MPI # put this first!
using StochJuMP, JuMP

MPI.init()

numScens = 2
m = StochasticModel(numScens)

@defVar(m, x >= 1)
@defVar(m, y <= 2)

@addConstraint(m,  x + y == 1)
@addConstraint(m, -x + y <= 1)
setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y)

comm = MPI.COMM_WORLD
size = MPI.size(comm)
rank = MPI.rank(comm)
scenPerRank = iceil(numScens/size)
proc_idx_set = rank*scenPerRank + (1:scenPerRank)
if endof(proc_idx_set) > numScens # handle case where numScens is not a multiple of size
    proc_idx_set = (rank*scenPerRank+1):numScens
end

println("rank = $rank, size = $size")
println("scenPerRank  = $scenPerRank")
println("proc_idx_set = $proc_idx_set")

rhs  = [5,4]
coef = [2,3]
qc   = [1,0.5]
ac   = [1,0.75]

for i in proc_idx_set
    bl = StochasticBlock(m)
    @defVar(bl, 0 <= w <= 1)
    @addConstraint(bl, coef[i]w - x - y <= rhs[i])
    @addConstraint(bl, coef[i]w + x     == rhs[i])
    setObjective(bl, :Min, qc[i]*w*w + ac[i]*w)
end

StochJuMP.pips_solve(m)
