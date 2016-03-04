# StochJuMP.jl
The StochJuMP.jl package provides a scalable algebraic modeling framework for stochastic programming problems in Julia. StochJuMP.jl is an extension of the [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl) package, which is as fast as [AMPL](http://ampl.com) and faster than any other modeling tools such as [GAMS](http://www.gams.com) and [Pyomo](http://www.pyomo.org) (see [this](http://arxiv.org/pdf/1312.1431.pdf)).

An example of the StochJuMP.jl package reads:
```julia
using StochJuMP

numScen = 2
m = StochasticModel(numScen)

@defVar(m, 0 <= x <= 1)
@defVar(m, 0 <= y <= 1)

@addConstraint(m, x + y == 1)
setObjective(m, :Min, x*x + y)

for i in 1:numScen
    bl = StochasticModel(parent=m)
    @defVar(bl, w >= 0)
    @addConstraint(bl, w - x - y <= 1)
    @setObjective(bl, Min, w*w + w)
end
```

## Solvers for StochJuMP
The StochJuMP model can be solved by either PIPS or DSP. [PIPS](http://git.mcs.anl.gov/PIPS.git/) is an open-source parallel interior point solver for stochastic convex and nonconvex continuous programs. [DSP](https://github.com/kibaekkim/DSP) is a open-source package of the parallel decomposition methdos for stochastic mixed-integer programs. The Julia interface for PIPS and DSP are also available in [PIPS.jl](https://github.com/kibaekkim/PIPS.jl) and [DSPsolver.jl](https://github.com/kibaekkim/DSPsolver.jl), respectively.
