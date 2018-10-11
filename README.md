# StructJuMP
The StructJuMP package provides a parallel algebraic modeling framework for block structured optimization models in Julia. StructJuMP, originally known as StochJuMP, is tailored to two-stage stochastic optimization problems and uses MPI to enable a parallel, distributed memory instantiation of the problem. StructJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).

[![Build Status](https://travis-ci.org/StructJuMP/StructJuMP.jl.svg?branch=master)](https://travis-ci.org/StructJuMP/StructJuMP.jl)

## Installation

The most recently tagged version of this package works with
most recently tagged version of [JuMP](https://github.com/JuliaOpt/JuMP.jl),
i.e. JuMP v0.18.x, to try this package with JuMP v0.18, simply do
```julia
] add StructJuMP
```
The development of this package for JuMP v0.18.x is done in the `release-0.1`
branch.
The `master` branch of this package works with the JuMP v0.19-alpha. To try this
package with JuMP v0.19-alpha, do:
```julia
] add JuMP#v0.19-alpha
] add StructJuMP#master
```

```julia
Pkg.add("StructJuMP")
```

## An example
```julia
using StructJuMP

numScen = 2
m = StructuredModel(num_scenarios=numScen)
@variable(m, x[1:2])
@NLconstraint(m, x[1] + x[2] == 100)
@NLobjective(m, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for i in 1:numScen
    bl = StructuredModel(parent=m, id=i)
    @variable(bl, y[1:2])
    @NLconstraint(bl, x[1] + y[1] + y[2] ≥  0)
    @NLconstraint(bl, x[2] + y[1] + y[2] ≤ 50)
    @NLobjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end
```
The above example builds a two level structured model `m` with 2 scenarios.

## Available Solvers for StructJuMP

### Nonlinear Solvers
Problems modeled in StructJuMP models can be solved in parallel using the [PIPS-NLP](https://github.com/Argonne-National-Laboratory/PIPS) parallel optimization solver. In addition, StructJuMP models can be solved (in serial only) using [Ipopt](https://projects.coin-or.org/Ipopt). The SturctJuMP models interface with the solvers via [StructJuMPSolverInterface.jl](https://github.com/Argonne-National-Laboratory/StructJuMPSolverInterface.jl).

### Mixed-Integer Solvers
[DSP](https://github.com/Argonne-National-Laboratory/DSP.git) can read models from StructJuMP via [DSPsolver.jl](https://github.com/Argonne-National-Laboratory/DSPsolver.jl.git). In particular, ``DSP`` can solver problems with integer variables in parallel.

### Stochastic Dual Dynamic Programming (SDDP)
[StructDualDynProg](https://github.com/blegat/StructDualDynProg.jl) can run the SDDP algorithm on multi-stage models from StructJuMP.

## Acknowledgements
StructJuMP has been developed under the financial support of Department of Energy (DOE), Office of Advanced Scientific Computing Research, Office of Electricity Delivery and Energy Reliability, and Grid Modernization Laboratory Consortium (GMLC) (PIs: Cosmin G. Petra, Lawrence Livermore National Laboratory and Mihai Anitescu, Argonne National Laboratory).
