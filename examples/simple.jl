using StructJuMP

m = StructuredModel(num_scenarios=p)

n = 10
@variable(m, x[1:n] >= 0)
@variable(m, y[1:n] >= 0)
for i in 1:n
    @constraint(m, x[i] - y[i] <= i)
end

p = 3
q = 2
for s = 1:p
    bl = StructuredModel(parent=m, num_scenarios=2, id=s)
    @variable(bl, z[1:n] >= 0)
    @constraint(bl, x + sum(z[i] for i=1:n) .== 1)
    for t = 1:q
        bll = StructuredModel(parent=bl, id=t)
        @variable(bll, w[1:n] >= 0)
        @constraint(bll, sum(w[i] for i=1:n if iseven(i)) == 1)
        par = parent(bll)
        @constraint(bll, sum(z[i] for i=1:n) - sum(w[i] for i=1:n) >= 0)
        @constraint(bll, sum(z[i] for i=1:n) >= 0)
        @constraint(bll, 0 >= sum(w[i] for i=1:n))
    end
end

@constraint(m, sum(z for bl in children(m), z=getVar(bl,:z)) == 1)


#sum(indx for cond(indx)) (models[indx].y[1] + x[1]) = 0
