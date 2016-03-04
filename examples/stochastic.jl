using StochJuMP

m = StochasticModel()

@defStochasticVar(m, xh[ASSETS] >= 0)
@defStochasticVar(m, xb[ASSETS] >= 0)
@defStochasticVar(m, xs[ASSETS] >= 0)

# Stage 0
sub = StochasticModel(parent=m)
@addConstraint(sub, (1+Ct)*sum{val[j]*xb[j], j = ASSETS} <= BUDGET)

for stage = 1:T
    sub = StochasticModel(parent=m) # return a ConstraintRef to a Stage?
    for j in ASSETS
        @addConstraint(sub, xh[j] == (1+Ret[j])*xh[getParent(sub),j] + xb[j] - xs[j])
        @addConstraint(sub, (1-Ct)*sum{val[j]*xs[j], j = ASSETS} == LIABILITY + (1+Ct)*sum{val[j]*xb[j], j = ASSETS})
    end
end

sub = StochasticModel(parent=m)
@defVar(sub, wealth)
@addConstraint(sub, wealth == sum{val[j]*xh[j], j=ASSETS})
@addConstraint(sub, @Exp(wealth) == mu)
setObjective(sub, mu - Rho*((wealth*wealth)- mu*mu))

# Exp(x) = sum{ndinNODESET:ndincurrentstage}Prob[nd] ∗x[nd]
