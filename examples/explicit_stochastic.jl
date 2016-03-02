using StochJuMP

######
# DATA
######
NUMSTAGES = 3
TESTTIME = GEN = WIND = BUS = LIN = GEN = LOAD = 1:2
NODES = HORIZON = [1]
FIXEDGENCOST = THETASCALE = COSTSCALE = 1
Pmax = zeros(2,2)
np_cap = wind_share = X = V = max_ur = Pgen_init = gen_cost = ones(2)
rec_bus = snd_bus = bus_gen = bus_load = bus_wind = [1:2]
wind_total = [1 1]
node = NODES[1]
Pload = ones(2,2)
ref_bus = 1

m = StochasticModel()

@defStochasticVar(m, dummyV >= 0)

s0 = StochasticModel(parent=m)
@defStochasticVar(s0, 0 <= y[TESTTIME,GEN] <= 1)
@setObjective(s0, Min, FIXEDGENCOST*sum{y[t,j], t=TESTTIME,j=GEN})

for it in 1:NUMSTAGES
    st = StochasticModel(parent=s0)
    @defStochasticVar(st, 0 <= Pgen[TESTTIME,j=GEN] <= np_cap[j])
    @defStochasticVar(st, 0 <= Pwind[t=TESTTIME,j=WIND] <= wind_total[node,t]*wind_share[j])
    @defStochasticVar(st, -THETASCALE*π/2 <= theta[TESTTIME,j=BUS] <= THETASCALE*π/2)
    @defStochasticVar(st, -Pmax[i] <= P[TESTTIME,i=LIN] <= Pmax[i])
    for t in TESTTIME
        for j in BUS
            @addConstraint(st, ( sum{P[t,i], i=LIN; j==rec_bus[i]}
                                -sum{P[t,i], i=LIN; j==snd_bus[i]}
                                +sum{Pgen[t,i], i=GEN; j==bus_gen[i]}
                                -sum{Pload[t,i], i=LOAD; j==bus_load[i]}
                                +sum{wind_total[node,t]*wind_share[i]-Pwind[t,i], i=WIND; j==bus_wind[i]}
                                == 0))
        end
        @addConstraint(st, theta[t,ref_bus] == 0)
        for i in LIN
            @addConstraint(st, X[i]*P[t,i] - V[i]^2*(theta[t,snd_bus[i]]-theta[t,rec_bus[i]])/THETASCALE == 0)
        end
    end
    # ramp constraints
    for t in HORIZON
        for i in GEN
            @addConstraint(st, -max_ur[i] <= Pgen[t,i] - Pgen[t+1,i] <= max_ur[i]) # range constraint, beware...
        end
    end

    for i in GEN
        @addConstraint(st, -max_ur[i] <= Pgen[1,i] - Pgen_init[i] <= max_ur[i])
    end

    # linking constraints
    # gen[s,t,j] <= np_cap[j]*y[t,j]
    for t in TESTTIME
        for j in GEN
            @addConstraint(st, -np_cap[j] <= Pgen[t,j] - np_cap[j]*y[t,j] <= 0) # should wrap the call to y in an ancestor-type thing
        end
    end

    # need Exp here for earlier sml versions
    @setObjective(st, Min, sum{COSTSCALE*gen_cost[i]*Pgen[t,i], t=TESTTIME, i=GEN})

end
