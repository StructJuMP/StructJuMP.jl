using StructJuMP

######
# DATA
######
NUMSCEN = 3
TESTTIME = GEN = WIND = BUS = LIN = GEN = LOAD = 1:2
NODES = HORIZON = [1]
FIXEDGENCOST = THETASCALE = COSTSCALE = 1
Pmax = zeros(2,2)
np_cap = wind_share = X = V = max_ur = Pgen_init = gen_cost = ones(2)
rec_bus = snd_bus = bus_gen = bus_load = bus_wind = 1:2
wind_total = [1 1]
node = NODES[1]
Pload = ones(2,2)
ref_bus = 1

m = StructuredModel()

@variable(m, dummyV >= 0)

s0 = StructuredModel(parent=m, id=1)
@variable(s0, 0 <= y[TESTTIME,GEN] <= 1)
@objective(s0, :Min, FIXEDGENCOST*sum(y[t,j] for t=TESTTIME,j=GEN))

for it=1:NUMSCEN
    st = StructuredModel(parent=s0, id=1+it)
    @variable(st, 0 <= Pgen[TESTTIME,j=GEN] <= np_cap[j])
    @variable(st, 0 <= Pwind[t=TESTTIME,j=WIND] <= wind_total[node,t]*wind_share[j])
    @variable(st, -THETASCALE*π/2 <= theta[TESTTIME,j=BUS] <= THETASCALE*π/2)
    @variable(st, -Pmax[i] <= P[TESTTIME,i=LIN] <= Pmax[i])
    for t in TESTTIME
        for j in BUS
            @constraint(st, ( sum(P[t,i] for i=LIN if j==rec_bus[i])
                                -sum(P[t,i] for i=LIN if j==snd_bus[i])
                                +sum(Pgen[t,i] for i=GEN if j==bus_gen[i])
                                -sum(Pload[t,i] for i=LOAD if j==bus_load[i])
                                +sum(wind_total[node,t]*wind_share[i]-Pwind[t,i] for i=WIND if j==bus_wind[i])
                                == 0))
        end
        @constraint(st, theta[t,ref_bus] == 0)
        for i in LIN
            @constraint(st, X[i]*P[t,i] - V[i]^2*(theta[t,snd_bus[i]]-theta[t,rec_bus[i]])/THETASCALE == 0)
        end
    end
    # ramp constraints
    for t in HORIZON
        for i in GEN
            @constraint(st, -max_ur[i] <= Pgen[t,i] - Pgen[t+1,i] <= max_ur[i]) # range constraint, beware...
        end
    end

    for i in GEN
        @constraint(st, -max_ur[i] <= Pgen[1,i] - Pgen_init[i] <= max_ur[i])
    end

    # linking constraints
    # gen[s,t,j] <= np_cap[j]*y[t,j]
    for t in TESTTIME
        for j in GEN
            @constraint(st, -np_cap[j] <= Pgen[t,j] - np_cap[j]*y[t,j] <= 0) # should wrap the call to y in an ancestor-type thing
        end
    end

    # need Exp here for earlier sml versions
    @objective(st, :Min, sum(COSTSCALE*gen_cost[i]*Pgen[t,i] for t=TESTTIME, i=GEN))

end
