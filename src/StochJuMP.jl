module StochJuMP

import MPI
using JuMP # To reexport, should be using (not import)
import MathProgBase
import MathProgBase.MathProgSolverInterface

export StochasticModel, getStochastic, getparent, getchildren, getProcIdxSet,
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


# --------------
# StochasticData
# --------------

# Teyp Define
type StochasticData
    probability::Vector{Number}
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
end

# Constructor with no argument
StochasticData() = StochasticData(Number[], JuMP.Model[], nothing, 0)

# Constructor without specifyng probabilities
StochasticData(children, parent, nscen) = StochasticData(Number[], children, parent, nscen)


# ---------------
# StochasticModel
# ---------------

# Constructor with the number of scenarios
function StochasticModel(numScen::Int)
    # MPI.init()
    m = JuMP.Model()
    m.ext[:Stochastic] = StochasticData(JuMP.Model[],nothing,numScen)
    return m
end

# Constructor with children and partent models
StochasticModel(children,parent) = StochasticModel(children,parent,0)
function StochasticModel(children, parent, nscen)
    m = JuMP.Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(children, parent, nscen)
    return m
end


# -------------
# Get functions
# -------------

function getStochastic(m::JuMP.Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

getparent(m::JuMP.Model)      = getStochastic(m).parent
getchildren(m::JuMP.Model)    = getStochastic(m).children
getprobability(m::JuMP.Model) = getStochastic(m).probability
num_scenarios(m::JuMP.Model)  = getStochastic(m).num_scen

function getProcIdxSet(m::JuMP.Model)
    numScens = num_scenarios(m)
    comm = MPI.COMM_WORLD
    mysize = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    # Why don't we just take a round-and-robin?
    proc_idx_set = Int[];
    for s = myrank:mysize:(numScens-1)
        push!(proc_idx_set, s+1);
    end
    return proc_idx_set;
end

# ---------------
# StochasticBlock
# ---------------

# Constructor without probability
StochasticBlock(m::JuMP.Model) = StochasticBlock(m::JuMP.Model, 1.0 / num_scenarios(m))

# Construcor with probability
function StochasticBlock(m::JuMP.Model, probability::Number)
    stoch = getStochastic(m)
    #stoch.parent = m
    ch = StochasticModel(JuMP.Model[], m)
    push!(stoch.children, ch)
    push!(stoch.probability, probability)
    return ch
end

macro second_stage(m,ind,code)
    return quote
        proc_idx_set = getProcIdxSet($(esc(m)))
        for $(esc(ind)) in proc_idx_set
            $(esc(code))
        end
    end
end

end
