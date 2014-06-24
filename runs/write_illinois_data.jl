using DataFrames

fp = open("$(ENV["HOME"])/.julia/v0.3/StochJuMP/runs/data.csv", "w")

NS_max = 256

df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Lines_data.tab", separator='\t', skipstart=1)
LIN          = df[:LIN]
snd_bus      = Dict(LIN,df[:snd_bus])
rec_bus      = Dict(LIN,df[:rec_bus])
Pmax         = Dict(LIN,df[:Pmax])
snd_bus_arr = Array(Int, length(LIN))
rec_bus_arr = Array(Int, length(LIN))
Pmax_arr    = Array(Float64, length(LIN))
for (it,lin) in enumerate(LIN)
    snd_bus_arr[it] = snd_bus[lin]
    rec_bus_arr[it] = rec_bus[lin]
    Pmax_arr[it]    = Pmax[lin]
end
for it in 1:(length(LIN)-1)
    print(fp, snd_bus_arr[it])
    print(fp, ", ")
end
println(fp, snd_bus_arr[end])
for it in 1:(length(LIN)-1)
    print(fp, rec_bus_arr[it])
    print(fp, ", ")
end
println(fp, rec_bus_arr[end])
for it in 1:(length(LIN)-1)
    print(fp, Pmax_arr[it])
    print(fp, ", ")
end
println(fp, Pmax_arr[end])

df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/bus_data.tab", separator='\t', skipstart=1)
BUS     = df[:BUS]
for i in 1:(length(BUS)-1)
    print(fp, BUS[i])
    print(fp, ", ")
end
println(fp, BUS[end])

# thermal generators
df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Gen_data_thermals.tab", separator='\t', skipstart=1)
GENTHE       = df[:GENTHE]
bus_genThe   = Dict(GENTHE,df[:bus_genThe])
np_capThe    = Dict(GENTHE,df[:np_capThe])
min_hrateThe = Dict(GENTHE,df[:min_hrateThe])
fuelThe      = Dict(GENTHE,df[:fuelThe])
bus_gen_arr = Array(Int, length(GENTHE))
np_cap_arr  = Array(Float64, length(GENTHE))
for (it,lin) in enumerate(GENTHE)
    bus_gen_arr[it] = bus_genThe[lin]
    np_cap_arr[it]  = np_capThe[lin]
end
for it in 1:(length(GENTHE)-1)
    print(fp, bus_gen_arr[it])
    print(fp, ", ")
end
println(fp, bus_gen_arr[end])
for it in 1:(length(GENTHE)-1)
    print(fp, np_cap_arr[it])
    print(fp, ", ")
end
println(fp, np_cap_arr[end])

# wind generators
df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Gen_data_wind.tab", separator='\t', skipstart=1)
GENWIN     = df[:GENWIN]
bus_genWin = Dict(GENWIN,df[:bus_genWin])
np_capWin  = Dict(GENWIN,df[:np_capWin])
fuelWin    = Dict(GENWIN,df[:fuelWin])
bus_gen_arr = Array(Int, length(GENWIN))
np_cap_arr  = Array(Float64, length(GENWIN))
for (it,lin) in enumerate(GENWIN)
    bus_gen_arr[it] = bus_genWin[lin]
    np_cap_arr[it]  = np_capWin[lin]
end
for it in 1:(length(GENWIN)-1)
    print(fp, bus_gen_arr[it])
    print(fp, ", ")
end
println(fp, bus_gen_arr[end])
for it in 1:(length(GENWIN)-1)
    print(fp, np_cap_arr[it])
    print(fp, ", ")
end
println(fp, np_cap_arr[end])

# fuels
df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/fuel_data_distinctPrices.tab", separator='\t', skipstart=1)
FUEL        = df[:FUEL]
HV          = Dict(FUEL,df[:HV])
Unitprice   = Dict(FUEL,df[:Unitprice])

# loads
df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/load_load.tab", separator='\t', skipstart=1)
LOAD     = df[:LOAD]
bus_load = Dict(LOAD,df[:bus_load])
bus_load_arr = Array(Int, length(LOAD))
for (it,lin) in enumerate(LOAD)
    bus_load_arr[it] = bus_load[lin]
end
for it in 1:(length(LOAD)-1)
    print(fp, bus_load_arr[it])
    print(fp, ", ")
end
println(fp, bus_load_arr[end])

df = readdlm("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/Loads.dat", '\t')
loads = df[3,:]
for i in 1:length(loads)
     loads[i] = (loads[i] > 1000 ? 1000 : 1.2*loads[i])
end
for i in 1:(length(loads)-1)
    print(fp, loads[i])
    print(fp, ", ")
end
println(fp, loads[end])

gen_cost_the = Array(Float64, length(GENTHE))
for (it,i) in enumerate(GENTHE)
     gen_cost_the[it] = 1e-3*min_hrateThe[i] / HV[fuelThe[i]] * Unitprice[fuelThe[i]]
end
for i in 1:(length(GENTHE)-1)
    print(fp, gen_cost_the[i])
    print(fp, ", ")
end
println(fp, gen_cost_the[end])

gen_cost_win = Array(Float64, length(GENWIN))
for (it,i) in enumerate(GENWIN)
     gen_cost_win[it] = HV[fuelWin[i]]*Unitprice[fuelWin[i]]
end
for i in 1:(length(GENWIN)-1)
    print(fp, gen_cost_win[i])
    print(fp, ", ")
end
println(fp, gen_cost_win[end])

df = readtable("$(ENV["HOME"])/.julia/v0.3/StochJuMP/examples/Illinois/IIDmean_2006_06_04_0_0.dat", header=false)
windPower = Array(Float64, length(df[:x1]))
for (it, val) in enumerate(df[:x1])
    windPower[it] = val
end
for i in 1:(length(windPower)-1)
    print(fp, windPower[i])
    print(fp, ", ")
end
println(fp, windPower[end])
for s in 2:NS_max
    for gw in 1:length(GENWIN)
        val = windPower[gw] + 0.25windPower[gw]*randn()
        val = min(10*(1+(exp(2*0.7*1.2*val-4)-1)/(exp(2*0.7*1.2*val-4)+1)),np_cap_arr[gw])
        print(fp, "$val, ")
    end
    println(fp, )
end
close(fp)
