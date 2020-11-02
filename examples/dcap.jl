using StructJuMP
using Random

function DCAP(nR::Int, nN::Int, nT::Int, nS::Int, seed::Int=1)::StructuredModel

    # set random seed (default=1)
    Random.seed!(seed)

    # generate & store instance data
    ## sets
    R = 1:nR
    N = 1:nN
    T = 1:nT
    S = 1:nS

    ## parameters
    a = rand(nR, nT) * 5 .+ 5
    b = rand(nR, nT) * 40 .+ 10
    c = rand(nR, nN, nT, nS) * 5 .+ 5
    c0 = rand(nN, nT, nS) * 500 .+ 500
    d = rand(nN, nT, nS) .+ 0.5
    Pr = ones(nS)/nS

    # construct JuMP.Model
    model = StructuredModel(num_scenarios = nS)

    ## 1st stage
    @variable(model, x[i=R,t=T] >= 0)
    @variable(model, u[i=R,t=T], Bin)
    @objective(model, Min, sum(a[i,t]*x[i,t] + b[i,t]*u[i,t] for i in R for t in T))
    @constraint(model, [i=R,t=T], x[i,t] - u[i,t] <= 0)

    ## 2nd stage
    for s in S
        sb = StructuredModel(parent=model, id = s, prob = Pr[s])
        # @variable(sb, y[i=R, j=N, t=T], Bin)
        @variable(sb, 0 <= y[i=R, j=N, t=T] <= 1)
        #@variable(sb, z[j=N,t=T] >= 0) # originally implemented variable (continuous)
        # @variable(sb, z[j=N,t=T], Bin)  # modify as SIPLIB 1.0
        @variable(sb, 0 <= z[j=N,t=T] <= 1)
        @objective(sb, Min, sum(c[i,j,t,s]*y[i,j,t] for i in R for j in N for t in T) + sum(c0[j,t,s]*z[j,t] for j in N for t in T))
        @constraint(sb, [i=R, t=T], -sum(x[i,tau] for tau in 1:t) + sum(d[j,t,s]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N, t=T], sum(y[i,j,t] for i in R) + z[j,t] == 1)
    end

    return model
end