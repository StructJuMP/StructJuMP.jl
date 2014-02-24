m = StochasticModel()

@defVar(m, xh[ASSETS] >= 0)
@defVar(m, xb[ASSETS] >= 0)
@defVar(m, xs[ASSETS] >= 0)

# Stage 0
sub = Stage(m, 0)
@addConstraint(sub, (1+Ct)*sum{val[j]*xb[j], j = ASSETS} <= BUDGET)

for stage = 1:T
    sub = Stage(m,stage) # return a ConstraintRef to a Stage?
    for j in ASSETS
        @addConstraint(sub, xh[j] == (1+Ret[j])*sub.ancestor[1].xh[j] + xb[j] - xs[j])
        @addConstraint(sub, (1-Ct)*sum{val[j]*xs[j], j = ASSETS} == LIABILITY + (1+Ct)*sum{val[j]*xb[j], j = ASSETS})
    end
end

sub = Stage(m,T)
@defVar(sub, wealth)
@addConstraint(sub, wealth == sum{val[j]*xh[j], j=ASSETS})
@addConstraint(sub, @Exp(wealth) == mu)
setObjective(sub, mu - Rho*((wealth*wealth)- mu*mu))

# Exp(x) = sum{ndinNODESET:ndincurrentstage}Prob[nd] ∗x[nd]
