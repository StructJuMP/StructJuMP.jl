using StochJuMP

m = StochasticModel()

n = 10
@defVar(m, x[1:n] >= 0)
@defVar(m, y[1:n] >= 0)
for i in 1:n
    @addConstraint(m, x[i] - y[i] <= i)
end

p = 3
q = 2
for s = 1:p
    bl = StochasticBlock(m, "level1-$s")
    @defVar(bl, z[1:n] >= 0)
    @addConstraint(bl, sum{z[i], i=1:n} == 1)
    for t = 1:q
        bll = StochasticBlock(bl, "level2-$s,$t")
        @defVar(bll, w[1:n] >= 0)
        @addConstraint(bll, sum{w[i], i=1:n; iseven(i)} == 1)
        par = parent(bll)
        @addConstraint(bll, sum{z[i], i=1:n} - sum{w[i], i=1:n} >= 0)
        @addConstraint(bll, sum{z[i], i=1:n} >= 0)
        @addConstraint(bll, 0 >= sum{w[i], i=1:n})
    end
end
