# StructJuMP
The StructJuMP package provides a parallel algebraic modeling framework for block structured optimization models in Julia. StructJuMP, originally known as StochJuMP, is tailored to two-stage stochastic optimization problems and uses MPI to enable a parallel, distributed memory instantiation of the problem. StructJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).


##Nonlinear solvers

Problems modeled in StructJuMP models can be solved in parallel using the [PIPS-NLP](https://github.com/Argonne-National-Laboratory/PIPS) parallel optimization solver. In addition, StructJuMP models can be solved (in serial only) using [Ipopt](https://projects.coin-or.org/Ipopt). The solver interface and glue code  for PIPS-NLP and Ipopt are located at [StructJuMPSolverInterface](https://github.com/Argonne-National-Laboratory/StructJuMPSolverInterface.jl). 


##Installation

Prior its installation, StructJuMP requires both PIPS-NLP and Ipopt solvers to be installed. 

1. Install and build the PIPS-NLP solver using the installation instructions found [here](https://github.com/Argonne-National-Laboratory/PIPS). Then set the environment variables `PIPS_NLP_PAR_SHARED_LIB` and `PIPS_NLP_SHARED_LIB` to the location of the PIPS-NLP parallel and serial shared (.so|.dylib) libraries. Note that these shared libraries are named `libparpipsnlp.so|.dylib` and `libpipsnlp.so|.dylib`.

2. Install the Ipopt solver, if not already installed:

 `julia> Pkg.add("Ipopt")`

3. Install StructJuMP (this package):

 `julia> Pkg.clone("https://github.com/joehuchette/StructJuMP.jl")

4. Install StructJuMPSolverInterface:
 
 `julia> Pkg.clone("https://github.com/Argonne-National-Laboratory/StructJuMPSolverInterface.jl")`


##An example
```julia
using StructJuMP, JuMP

numScen = 2
m = StructuredModel(num_scenarios=numScen)
@variable(m, x[1:2])
@NLconstraint(m, x[1] + x[2] == 100)
@NLobjective(m, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for i in 1:numScen
    bl = StructuredModel(parent=m, id=i)
    @variable(bl, y[1:2])
    @NLconstraint(bl, x[1] + y[1]+y[2] ≥  0)
    @NLconstraint(bl, x[2] + y[1]+y[2] ≤ 50)
    @NLobjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end
```
The above example builds a two level structured model `m` with 2 scenarios. The model can be solved by calling solve function with a solver that implements the StructJuMPSolverInterface. 
```julia
solve(m, solver="PipsNlp") #solving using parallel PIPS-NLP solver
solve(m, solver="PipsNlpSerial") #solving using serial PIPS-NLP solver, mostly for debug and maintenance purpose
solve(m, solver="Ipopt") #solving using (serial) Ipopt solver
```

## Known Limitation
* If a constraint declared at the sub-problem uses variable from the parent level, it has to be add using @NLconstraint (instead of @constraint). 


## Acknowledgements
StructJuMP has been developed under the financial support of Department of Energy (DOE), Office of Advanced Scientific Computing Research, Office of Electricity Delivery and Energy Reliability, and Grid Modernization Laboratory Consortium (GMLC) (PIs: Cosmin G. Petra and Mihai Anitescu, Argonne National Laboratory).
