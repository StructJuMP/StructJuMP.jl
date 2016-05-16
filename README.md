# StructJuMP.jl
The StructJuMP.jl package provides a scalable algebraic modeling framework for block structured optimization models in Julia. StructJuMP was originally known as StochJuMP, and tailored specifically towards two-stage stochastic optimization problems. StructJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).

An example of the StructJuMP.jl package reads:
```julia
using StructJuMP, JuMP

scen = 20
m = StructuredModel(num_scenarios=scen)
@defVar(m, x[1:2])
@addNLConstraint(m, x[1] + x[2] == 100)
@setNLObjective(m, Min, x[1]^2 + x[2]^2 + x[1]*x[2])

for i in 1:scen
    bl = StructuredModel(parent=m)
    @defVar(bl, y[1:2])
    @addNLConstraint(bl, x[1] + y[1]+y[2] ≥  0)
    @addNLConstraint(bl, x[2] + y[1]+y[2] ≤ 50)
    @setNLObjective(bl, Min, y[1]^2 + y[2]^2 + y[1]*y[2])
end
```

## Solvers for StructJuMP
### Structured solver
The StructJuMP model can be solved by either PIPS or DSP. [PIPS](https://github.com/Argonne-National-Laboratory/PIPS/) is an open-source parallel interior point solver for stochastic convex and nonconvex continuous programs. [DSP](https://github.com/kibaekkim/DSP) is a open-source package of the parallel decomposition methods for stochastic mixed-integer programs. The Julia interface for PIPS and DSP are also available in [PIPS.jl](https://github.com/kibaekkim/PIPS.jl) and [DSPsolver.jl](https://github.com/kibaekkim/DSPsolver.jl), respectively.
Example solver interface implementations for PIPS are provided in [PIPS Structure Interface](https://github.com/fqiang/StructJuMP.jl/blob/master/src/pips_structure_interface.jl) and [PIPS Serial Interface](https://github.com/fqiang/StructJuMP.jl/blob/master/src/serial_pipsnlp_interface.jl).

### Nonstructured solver
The StructJuMP model can also be solve by Ipopt. [Ipopt](https://projects.coin-or.org/Ipopt) is an Interior point optimization solver for finding local solution of large-scale nonlinear optimization problems. In this case, modeller can still take advantage of StructJuMP to model the problem using the problem' structure, but to solve the problem using more matured (or robust) solvers. An example interface implementation for Ipopt is provided in [Ipopt Interface](https://github.com/fqiang/StructJuMP.jl/blob/master/src/ipopt_interface.jl).

## Known Limitation
* If the constraint declared at the sub-problem uses variable from the parent level, it has to be add using @addNLConstraint. 
* Variables declared in sub-problem has to be putted in front of constraints declarations. This is limitation only applys when using the provided Ipopt and PIPS solver interface.  
