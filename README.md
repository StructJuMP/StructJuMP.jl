# StructJuMP.jl
The StructJuMP.jl package provides a scalable algebraic modeling framework for block structured optimization models in Julia. StructJuMP was originally known as StochJuMP, and tailored specifically towards two-stage stochastic optimization problems. StructJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).

An example of the StructJuMP.jl package reads:
```julia
using StructJuMP, JuMP

numScen = 2
m = StructuredModel(num_scenarios=numScen)
@variable(m, x[1:2])
@NLconstraint(m, x[1] + x[2] == 100)
@NLobjective(m, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for i in 1:scen
    bl = StructuredModel(parent=m)
    @variable(bl, y[1:2])
    @NLconstraint(bl, x[1] + y[1]+y[2] ≥  0)
    @NLconstraint(bl, x[2] + y[1]+y[2] ≤ 50)
    @NLobjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end
```

## Solvers for StructJuMP
The StructJuMP model can be solved by either PIPS or DSP. [PIPS](https://github.com/Argonne-National-Laboratory/PIPS/) is an open-source parallel interior point solver for stochastic convex and nonconvex continuous programs. [DSP](https://github.com/kibaekkim/DSP) is a open-source package of the parallel decomposition methods for stochastic mixed-integer programs. The Julia interface for PIPS and DSP are also available in [PIPS.jl](https://github.com/kibaekkim/PIPS.jl) and [DSPsolver.jl](https://github.com/kibaekkim/DSPsolver.jl), respectively.

##Nonlinear solvers for StructJuMP
The nonlinear StructJuMP models can be solved by either parallel or serial optimizaton solvers, such as PIPS-NLP and Ipopt respectively. [PIPS-NLP](https://github.com/Argonne-National-Laboratory/PIPS) is an extension to PIPS that parallelizes the solving algorithm for nonlinear programming problems using *Schur complement* approach. [Ipopt](https://projects.coin-or.org/Ipopt) is an Interior point optimization solver for finding local solution of large-scale nonlinear optimization problems. In the later case, modeller can still take advantage of StructJuMP to model the problem using its structure, and to solve the problem using more matured (or robust) solvers. The corresponding solver interface implementations for PIPS-NLP and Ipopt are located at [StructJuMPSolverInterface](https://github.com/fqiang/StructJuMPSolverInterface.jl). 

## Known Limitation
* If a constraint declared at the sub-problem uses variable from the parent level, it has to be add using @addNLConstraint (instead of @addConstraint). 
