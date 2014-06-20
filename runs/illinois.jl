import MPI
using JuMP, StochJuMP, DataFrames, Distributions

MPI.init()
myrank = MPI.rank(MPI.COMM_WORLD)
mysize = MPI.size(MPI.COMM_WORLD)

numScens = convert(Int, ARGS[1])

function solve_illinois(NS::Int)

    tic()
    SCEN = 1:NS
    NODES = 1:NS
    # lines
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Lines_data.tab", separator='\t', skipstart=1)
    LIN          = df[:LIN]
    snd_bus      = Dict(LIN,df[:snd_bus])
    rec_bus      = Dict(LIN,df[:rec_bus])
    Pmax         = Dict(LIN,df[:Pmax])

    # buses
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/bus_data.tab", separator='\t', skipstart=1)
    BUS     = df[:BUS]

    # thermal generators
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Gen_data_thermals.tab", separator='\t', skipstart=1)
    GENTHE       = df[:GENTHE]
    bus_genThe   = Dict(GENTHE,df[:bus_genThe])
    np_capThe    = Dict(GENTHE,df[:np_capThe])
    min_hrateThe = Dict(GENTHE,df[:min_hrateThe])
    fuelThe      = Dict(GENTHE,df[:fuelThe])

    # wind generators
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Gen_data_wind.tab", separator='\t', skipstart=1)
    GENWIN     = df[:GENWIN]
    bus_genWin = Dict(GENWIN,df[:bus_genWin])
    np_capWin  = Dict(GENWIN,df[:np_capWin])
    fuelWin    = Dict(GENWIN,df[:fuelWin])

    # fuels
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/fuel_data_distinctPrices.tab", separator='\t', skipstart=1)
    FUEL        = df[:FUEL]
    HV          = Dict(FUEL,df[:HV])
    Unitprice   = Dict(FUEL,df[:Unitprice])

    # loads
    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/load_load.tab", separator='\t', skipstart=1)
    LOAD     = df[:LOAD]
    bus_load = Dict(LOAD,df[:bus_load])

    df = readdlm("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Loads.dat", '\t')
    loads = df[3,:]
    for i in 1:length(loads)
         loads[i] = (loads[i] > 1000 ? 1000 : 1.2*loads[i])
    end

    gen_cost_the = Dict{Int,Float64}()
    for i in GENTHE
         gen_cost_the[i] = 1e-3*min_hrateThe[i] / HV[fuelThe[i]] * Unitprice[fuelThe[i]]
    end

    gen_cost_win = Dict{Int,Float64}()
    for i in GENWIN
         gen_cost_win[i] = HV[fuelWin[i]]*Unitprice[fuelWin[i]]
    end

    lineCutoff = 1

    df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/IIDmean_2006_06_04_0_0.dat", header=false)

    # generate windpower data
    windPower = Array(Dict{Int,Float64}, NS)
    windPower[1] = Dict(GENWIN,df[:x1])
    for s in 2:NS
        windPower[s] = Dict{Int,Float64}()
        for gw in GENWIN
            dist = Normal(windPower[1][gw],0.25*windPower[1][gw])
            windPower[s][gw] = rand(dist)
        end
    end
    for s in SCEN, gw in GENWIN
        windPower[s][gw] = min(10*(1+(exp(2*0.7*1.2*windPower[s][gw]-4)-1)/(exp(2*0.7*1.2*windPower[s][gw]-4)+1)),np_capWin[gw])
    end

    # model the thing
    m = StochasticModel(NS)

    # Stage 0
    @defVar(m, 0 <= Pgen_f[i=GENTHE] <= np_capThe[i])
    @defVar(m, 0 <= PgenWin_f[i=GENWIN] <= np_capWin[i])
    @defVar(m, -lineCutoff*Pmax[i] <= P_f[i=LIN] <= lineCutoff*Pmax[i])

    # (forward) power flow equations
    @addConstraint(m, pfeq_f[j=BUS],
                   +sum{P_f[i], i=LIN; j==rec_bus[i]}
                   -sum{P_f[i], i=LIN; j==snd_bus[i]}
                   +sum{Pgen_f[i], i=GENTHE; j==bus_genThe[i]}
                   +sum{PgenWin_f[i], i=GENWIN; j==bus_genWin[i]}
                   -sum{loads[i], i=LOAD; j==bus_load[i]} >= 0)

    #for node in NODES
    @second_stage m node begin
        bl = StochasticBlock(m)
        # variables
        @defVar(bl, 0 <= Pgen[i=GENTHE] <= np_capThe[i])
        @defVar(bl, 0 <= PgenWin[i=GENWIN] <= windPower[node][i])
        @defVar(bl, -lineCutoff*Pmax[i] <= P[i=LIN] <= lineCutoff*Pmax[i])

        # @addConstraint(bl, rampUp[g=GENTHE],
        #                Pgen[g] - Pgen_f[g] <=  np_capThe[g]/10)

        # @addConstraint(bl, randDown[g=GENTHE],
        #                Pgen[g] - Pgen_f[g] >= -np_capThe[g]/10)

        @addConstraint(bl, rampUpDown[g=GENTHE],
                       -np_capThe[g]/10 <= Pgen[g] - Pgen_f[g] <=  np_capThe[g]/10)

        # (spot) power flow equations
        @addConstraint(bl, pfeq[j=BUS],
                       +sum{P[i]-P_f[i], i=LIN; j==rec_bus[i]}
                       -sum{P[i]-P_f[i], i=LIN; j==snd_bus[i]}
                       +sum{Pgen[i]-Pgen_f[i], i=GENTHE; j==bus_genThe[i]}
                       +sum{PgenWin[i]-PgenWin_f[i], i=GENWIN; j==bus_genWin[i]} >= 0)

        @defVar(bl, t[GENTHE] >= 0)
        @addConstraint(bl, t_con1[g=GENTHE],
                        t[g] >= gen_cost_the[g]*Pgen_f[g] +
                        1.2*gen_cost_the[g]*(Pgen[g]-Pgen_f[g]))
        @addConstraint(bl, t_con2[g=GENTHE],
                        t[g] >= gen_cost_the[g]*Pgen_f[g])

        @defVar(bl, tw[GENWIN] >= 0)
        @addConstraint(bl, t_w_con1[g=GENWIN],
                        tw[g] >= gen_cost_win[g]*PgenWin_f[g] +
                        1.2*gen_cost_win[g]*(PgenWin[g]-PgenWin_f[g]))
        @addConstraint(bl, t_w_con2[g=GENWIN],
                        tw[g] >= gen_cost_win[g]*PgenWin_f[g])

        @setObjective(bl, Min, sum{ t[g], g=GENTHE} + sum{tw[g], g=GENWIN})
    end

    pips_time = StochJuMP.pips_solve(m)
    elapsed = toc()
    jump_time = elapsed - pips_time

    return jump_time, pips_time
end

# dummy call to compile everything
_,_ = solve_illinois(1)

jump_time, pips_time = solve_illinois(numScens)
if myrank == 0
    fp = open("$(ENV["HOME"])/.julia/v0.3/StochJuMP/runs/data/results.csv", "w")
    println(mysize, numScens, jump_time, pips_time)
    close(fp)
end

MPI.finalize()

