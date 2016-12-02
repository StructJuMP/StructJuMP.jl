using StructJuMP

m = StructuredModel()

@defStochasticVar(m, xh[ASSETS] >= 0)
@defStochasticVar(m, xb[ASSETS] >= 0)
@defStochasticVar(m, xs[ASSETS] >= 0)

# Stage 0
sub = StructuredModel(parent=m)
@constraint(sub, (1+Ct)*sum{val[j]*xb[j], j = ASSETS} <= BUDGET)

for stage = 1:T
    sub = StructuredModel(parent=m) # return a ConstraintRef to a Stage?
    for j in ASSETS
        @constraint(sub, xh[j] == (1+Ret[j])*xh[getParent(sub),j] + xb[j] - xs[j])
        @constraint(sub, (1-Ct)*sum{val[j]*xs[j], j = ASSETS} == LIABILITY + (1+Ct)*sum{val[j]*xb[j], j = ASSETS})
    end
end

sub = StructuredModel(parent=m)
@variable(sub, wealth)
@constraint(sub, wealth == sum{val[j]*xh[j], j=ASSETS})
@constraint(sub, @Exp(wealth) == mu)
JuMP.setobjective(sub, mu - Rho*((wealth*wealth)- mu*mu))

# Exp(x) = sum{ndinNODESET:ndincurrentstage}Prob[nd] ∗x[nd]
