using StructJuMP

m = StructuredModel()

n = 10
@variable(m, x[1:n] >= 0)
@variable(m, y[1:n] >= 0)
for i in 1:n
    @constraint(m, x[i] - y[i] <= i)
end

p = 3
q = 2
for s = 1:p
    bl = StructuredModel(parent=m)
    @variable(bl, z[1:n] >= 0)
    @constraint(bl, x + sum{z[i], i=1:n} == 1)
    for t = 1:q
        bll = StructuredModel(parent=bl)
        @variable(bll, w[1:n] >= 0)
        @constraint(bll, sum{w[i], i=1:n; iseven(i)} == 1)
        par = parent(bll)
        @constraint(bll, sum{z[i], i=1:n} - sum{w[i], i=1:n} >= 0)
        @constraint(bll, sum{z[i], i=1:n} >= 0)
        @constraint(bll, 0 >= sum{w[i], i=1:n})
    end
end

@constraint(m, sum{z, bl in children(m), z=getVar(bl,:z)} == 1)


sum{indx; cond(indx)} (models[indx].y[1] + x[1]) = 0
