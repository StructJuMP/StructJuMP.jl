import MPI
using JuMP, StochJuMP, DataFrames, Distributions

MPI.init()
comm = MPI.COMM_WORLD
myrank = MPI.rank(comm)
mysize = MPI.size(comm)

numScens = int(ARGS[1])
const root = 0

const LIN = 1:2522
const GENTHE = 1:193
const GENWIN = 1:870
const LOAD = 1:870

const lineCutoff = 1

function solve_illinois(NS::Int)
    tic()

    println("rank = $myrank, before data")
    if myrank == root
        fp = open("$(ENV["HOME"])/.julia/v0.3/StochJuMP/runs/data.csv", "r")
        line = chomp(readline(fp))
        snd_bus = int(split(line, ","))
        line = chomp(readline(fp))
        rec_bus = int(split(line, ","))
        line = chomp(readline(fp))
        Pmax = float(split(line, ","))
        line = chomp(readline(fp))
        bus_genThe = int(split(line, ","))
        line = chomp(readline(fp))
        np_capThe = float(split(line, ","))
        line = chomp(readline(fp))
        bus_genWin = int(split(line, ","))
        line = chomp(readline(fp))
        np_capWin = float(split(line, ","))
        line = chomp(readline(fp))
        bus_load = int(split(line, ","))
        line = chomp(readline(fp))
        loads = float(split(line, ","))
        line = chomp(readline(fp))
        gen_cost_the = float(split(line, ","))
        line = chomp(readline(fp))
        gen_cost_win = float(split(line, ","))
        windPower = Array(Float64, NS, length(GENWIN))
        line = chomp(readline(fp))
        windPower[1,:] = float(split(line, ","))
        for s in 1:NS
            line = chomp(readline(fp))
            windPower[s,:] = float(split(line, ","))
        end
        close(fp)
    else
        snd_bus      = Array(Int64,   length(LIN))
        rec_bus      = Array(Int64,   length(LIN))
        Pmax         = Array(Float64, length(LIN))
        bus_genThe   = Array(Int64,   length(GENTHE))
        np_capThe    = Array(Float64, length(GENTHE))
        bus_genWin   = Array(Int64,   length(GENWIN))
        np_capWin    = Array(Float64, length(GENWIN))
        bus_load     = Array(Int64,   length(LOAD))
        loads        = Array(Float64, length(LOAD))
        gen_cost_the = Array(Float64, length(GENTHE))
        gen_cost_win = Array(Float64, length(GENWIN))
        windPower    = Array(Float64, NS, length(GENWIN))
    end
    println("rank = $myrank, after data")
    MPI.Bcast!(snd_bus,      length(LIN),    root, comm)
    MPI.Bcast!(rec_bus,      length(LIN),    root, comm)
    MPI.Bcast!(Pmax,         length(LIN),    root, comm)
    MPI.Bcast!(bus_genThe,   length(GENTHE), root, comm)
    MPI.Bcast!(np_capThe,    length(GENTHE), root, comm)
    MPI.Bcast!(bus_genWin,   length(GENWIN), root, comm)
    MPI.Bcast!(np_capWin,    length(GENWIN), root, comm)
    MPI.Bcast!(bus_load,     length(LOAD),   root, comm)
    MPI.Bcast!(loads,        length(LOAD),   root, comm)
    MPI.Bcast!(gen_cost_the, length(GENTHE), root, comm)
    MPI.Bcast!(gen_cost_win, length(GENWIN), root, comm)
    for s in NS
        windPower[s,:] = MPI.Bcast!(windPower, length(GENWIN), root, comm)
    end
    println("rank = $myrank, after bcast")

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

    @second_stage m node begin
        bl = StochasticBlock(m)
        # variables
        @defVar(bl, 0 <= Pgen[i=GENTHE] <= np_capThe[i])
        @defVar(bl, 0 <= PgenWin[i=GENWIN] <= windPower[node][i])
        @defVar(bl, -lineCutoff*Pmax[i] <= P[i=LIN] <= lineCutoff*Pmax[i])

        @addConstraint(bl, rampUpDown[g=GENTHE],
                       -0.1np_capThe[g] <= Pgen[g] - Pgen_f[g] <=  0.1np_capThe[g])

        # (spot) power flow equations
        @addConstraint(bl, pfeq[j=BUS],
                       +sum{P[i]-P_f[i], i=LIN; j==rec_bus[i]}
                       -sum{P[i]-P_f[i], i=LIN; j==snd_bus[i]}
                       +sum{Pgen[i]-Pgen_f[i], i=GENTHE; j==bus_genThe[i]}
                       +sum{PgenWin[i]-PgenWin_f[i], i=GENWIN; j==bus_genWin[i]} >= 0)

        @defVar(bl, t[GENTHE] >= 0)
        @addConstraint(bl, t_con1[g=GENTHE],
                        t[g] >= gen_cost_the[g]*Pgen_f[g] +
                        1.2gen_cost_the[g]*(Pgen[g]-Pgen_f[g]))
        @addConstraint(bl, t_con2[g=GENTHE],
                        t[g] >= gen_cost_the[g]*Pgen_f[g])

        @defVar(bl, tw[GENWIN] >= 0)
        @addConstraint(bl, t_w_con1[g=GENWIN],
                        tw[g] >= gen_cost_win[g]*PgenWin_f[g] +
                        1.2gen_cost_win[g]*(PgenWin[g]-PgenWin_f[g]))
        @addConstraint(bl, t_w_con2[g=GENWIN],
                        tw[g] >= gen_cost_win[g]*PgenWin_f[g])

        @setObjective(bl, Min, sum{t[g], g=GENTHE} + sum{tw[g], g=GENWIN})
    end

    pips_time = StochJuMP.pips_solve(m)
    elapsed = toc()
    jump_time = elapsed - pips_time

    return jump_time, pips_time
end

# dummy call to compile everything
_,_ = solve_illinois(0)

jump_time, pips_time = solve_illinois(numScens)
if myrank == 0
    fp = open("$(ENV["HOME"])/.julia/v0.3/StochJuMP/runs/data/results.csv", "w")
    println(mysize, numScens, jump_time, pips_time)
    close(fp)
end

MPI.finalize()
