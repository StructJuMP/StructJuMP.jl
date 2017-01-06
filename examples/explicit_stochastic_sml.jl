using StructJuMP

m = StructuredModel()

@variable(m, xh0[ASSETS] >= 0)
@variable(m, xb0[ASSETS] >= 0)
@constraint(m, (1+Ct)*sum(val[j]*xb0[j] for j=ASSETS) <= BUDGET)

for n1 in NODES
    if PARENT[n1] == 0
        sub = StructuredModel(parent=m)
        @variable(sub, xh1[ASSETS] >= 0)
        @variable(sub, xb1[ASSETS] >= 0)
        @variable(sub, xs1[ASSETS] >= 0)
        for j in ASSETS
            @constraint(sub, xh1[j] == (1+Ret[n1,j])*xh0[j]+xb1[j]-xs1[j])
        end
        @constraint(sub, (1-Ct)*sum(val[j]*xs1[j] for j=ASSETS) == Liability1 + (1+Ct)*sum(val[j]*xb1[j]))
    end

    for n2 in NODES
        if PARENT[n2] == n1
            sub = StructuredModel(parent=m)
            @variable(sub, xh2[ASSETS] >= 0)
            @variable(sub, xb2[ASSETS] >= 0)
            @variable(sub, xs2[ASSETS] >= 0)
            for j in ASSETS
                @constraint(sub, xh2[j] == (1+Ret2[n2,j])*xh1[j] + xb2[j] - xs2[j])
            end
            @constraint(sub, (1-Ct)*sum(val[j]*xs2[j] for j=ASSETS) == Liability2 + (1+Ct)*sum(val[j]*xb2[j]))
            @variable(sub, wealth)
            @constraint(sub, wealth == Prob[n1]*Prob[n2]*sum(val[j]*xh2[j] for j=ASSETS))
            @objective(sub, mu - Rho*((wealth^2-mu^2))*Prob[n1]*Prob[n2])
        end
    end
end

@constraint(m, sum(Prob[n1]*Prob[n2]*wealth[(n1,n2)] for n1=NODES n2=NODES if PARENT[n1]==0 && PARENT[n2]==n1))

# things of note:
# Block(parent::Union(Model,Block), index) constructor
# should do name mangling so that variables get associated with their block
# add function block(parent, child1, child2, ...) that returns reference to subproblem. This should then be passable as index to subproblem variables
