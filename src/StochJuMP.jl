module StochJuMP

import MPI # think this needs to go first for MPI to work properly...

import JuMP

import MathProgBase
import MathProgBase.MathProgSolverInterface

export StochasticModel, getStochastic, getparent, getchildren, 
       num_scenarios, StochasticBlock, @second_stage

#############################################################################
# JuMP rexports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint, MultivarDict,
    ConstraintRef,
# Functions
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addSOS1, addSOS2, solve,
    getInternalModel, setPresolve, buildInternalModel,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    
# Macros and support functions
    @addConstraint, @addConstraints, @defVar, 
    @defConstrRef, @setObjective, addToExpression,
    @setNLObjective, @addNLConstraint, @gendict

type StochasticData
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
end

StochasticData() = StochasticData(JuMP.Model[],nothing,0)

function StochasticModel(numScen::Int)
    # MPI.init()
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
    return ch
end

macro second_stage(m,ind,code)
    return quote
        numScens = num_scenarios($(esc(m)))
        comm = MPI.COMM_WORLD
        size = MPI.size(comm)
        rank = MPI.rank(comm)
        scenPerRank = iceil(numScens/size)
        proc_idx_set = rank*scenPerRank + (1:scenPerRank)
        if endof(proc_idx_set) > numScens # handle case where numScens is not a multiple of size
            proc_idx_set = start(proc_idx_set):numScens
        end
        for $(esc(ind)) in proc_idx_set
            $(esc(code))
        end
    end
end

include("pips.jl")

end
