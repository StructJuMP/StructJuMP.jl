using StructJuMP, MPI

numScens = 2
m = StructuredModel(numScens)

@variable(m, x >= 1)
@variable(m, y <= 2)

@constraint(m,  x + y == 1)
@constraint(m, -x + y <= 1)
JuMP.setobjective(m, :Min, x*x + 0.5x*y + 0.25y*y - y)

rhs  = [5,4]
coef = [2,3]
qc   = [1,0.5]
ac   = [1,0.75]

@second_stage m scen begin
    bl = StructuredModel(parent=m)
    @variable(bl, 0 <= w <= 1)
    @constraint(bl, coef[scen]w - x - y <= rhs[scen])
    @constraint(bl, coef[scen]w + x     == rhs[scen])
    JuMP.setobjective(bl, :Min, qc[scen]*w*w + ac[scen]*w)
end
