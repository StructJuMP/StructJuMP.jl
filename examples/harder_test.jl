using StructJuMP, MPI

numScens = 2
m = StructuredModel(numScens)

@defVar(m, x >= 1)
@defVar(m, y <= 2)

@addConstraint(m,  x + y == 1)
@addConstraint(m, -x + y <= 1)
setObjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y)

rhs  = [5,4]
coef = [2,3]
qc   = [1,0.5]
ac   = [1,0.75]

@second_stage m scen begin
    bl = StructuredModel(parent=m)
    @defVar(bl, 0 <= w <= 1)
    @addConstraint(bl, coef[scen]w - x - y <= rhs[scen])
    @addConstraint(bl, coef[scen]w + x     == rhs[scen])
    setObjective(bl, :Min, qc[scen]*w*w + ac[scen]*w)
end
