# StructJuMP.jl
The StructJuMP.jl package provides a scalable algebraic modeling framework for block structured optimization models in Julia. StructJuMP was originally known as StochJuMP, and tailored specifically towards two-stage stochastic optimization problems. StructJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).

An example of the StructJuMP.jl package reads:
```julia
using StructJuMP

numScen = 2
m = StructuredModel(num_scenarios=numScen)

@defVar(m, 0 <= x <= 1)
@defVar(m, 0 <= y <= 1)

@addConstraint(m, x + y == 1)
setObjective(m, :Min, x*x + y)

for i in 1:numScen
    bl = StructuredModel(parent=m)
    @defVar(bl, w >= 0)
    @addConstraint(bl, w - x - y <= 1)
    @setObjective(bl, Min, w*w + w)
end
```

## Solvers for StructJuMP
The StructJuMP model can be solved by either PIPS or DSP. [PIPS](http://git.mcs.anl.gov/PIPS.git/) is an open-source parallel interior point solver for stochastic convex and nonconvex continuous programs. [DSP](https://github.com/kibaekkim/DSP) is a open-source package of the parallel decomposition methods for stochastic mixed-integer programs. The Julia interface for PIPS and DSP are also available in [PIPS.jl](https://github.com/kibaekkim/PIPS.jl) and [DSPsolver.jl](https://github.com/kibaekkim/DSPsolver.jl), respectively.
