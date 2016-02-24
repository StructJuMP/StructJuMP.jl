module StochJuMP

# import MPI
import JuMP # To reexport, should be using (not import)
import MathProgBase
import MathProgBase.MathProgSolverInterface
import ReverseDiffSparse

export StochasticModel, getStochastic, getparent, getchildren, getProcIdxSet,
       num_scenarios, StochasticBlock, @second_stage

# --------------
# StochasticData
# --------------

# Type Define
type StochasticData{T<:Number}
    probability::Vector{T}
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
    othervars::Vector{JuMP.Variable}
end


# Constructor with no argument
StochasticData() = StochasticData(Float64[], JuMP.Model[], nothing, 0, JuMP.Variable[])

# Constructor without specifyng probabilities
StochasticData(children, parent, nscen) = StochasticData(Float64[], children, parent, nscen, JuMP.Variable[])


# ---------------
# StochasticModel
# ---------------

# Constructor with the number of scenarios
function StochasticModel(;solver=JuMP.UnsetSolver(), num_scenarios::Int=0)
    m = JuMP.Model(solver=solver)
    m.ext[:Stochastic] = StochasticData(JuMP.Model[],nothing,num_scenarios)
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
        return m.ext[:Stochastic]::StochasticData
    else
        error("This functionality is only available for StochasticModels")
    end
end

getparent(m::JuMP.Model)      = getStochastic(m).parent
getchildren(m::JuMP.Model)    = getStochastic(m).children
getprobability(m::JuMP.Model) = getStochastic(m).probability
num_scenarios(m::JuMP.Model)  = getStochastic(m).num_scen

function getProcIdxSet(numScens::Integer)
    mysize = 1;
    myrank = 0;
    if isdefined(:MPI) == true && MPI.Initialized() == true
        comm = MPI.COMM_WORLD
        mysize = MPI.Comm_size(comm)
        myrank = MPI.Comm_rank(comm)
    end
    # Why don't we just take a round-and-robin?
    proc_idx_set = Int[];
    for s = myrank:mysize:(numScens-1)
        push!(proc_idx_set, s+1);
    end
    return proc_idx_set;
end

function getProcIdxSet(m::JuMP.Model)
    numScens = num_scenarios(m)
    return getProcIdxSet(numScens);
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

include("nlp.jl")

end
