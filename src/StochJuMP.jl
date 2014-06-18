module StochJuMP

import MPI # think this needs to go first for MPI to work properly...

import JuMP

import MathProgBase
import MathProgBase.MathProgSolverInterface

export StochasticModel, getStochastic, getparent, getchildren, num_scenarios, StochasticBlock

type StochasticData
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
end

StochasticData() = StochasticData(JuMP.Model[],nothing,0)

function StochasticModel(numScen::Int)
    m = JuMP.Model()
    m.ext[:Stochastic] = StochasticData(JuMP.Model[],nothing,numScen)
    return m
end

StochasticModel(children,parent) = StochasticModel(children,parent,0)
function StochasticModel(children, parent,nscen)
    m = JuMP.Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(children, parent, nscen)
    return m
end

function getStochastic(m::JuMP.Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

getparent(m::JuMP.Model)     = getStochastic(m).parent
getchildren(m::JuMP.Model)   = getStochastic(m).children
num_scenarios(m::JuMP.Model) = getStochastic(m).num_scen

function StochasticBlock(m::JuMP.Model)
    stoch = getStochastic(m)
    ch = StochasticModel(JuMP.Model[], m)
    push!(stoch.children, ch)
    stoch.num_scen += 1
    return ch
end

include("pips.jl")

end
