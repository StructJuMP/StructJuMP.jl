module StochJuMP

# import MPI
import JuMP # To reexport, should be using (not import)
import MathProgBase
import MathProgBase.MathProgSolverInterface
import ReverseDiffSparse

export StochasticModel, getStochastic, getparent, getchildren, getProcIdxSet,
       num_scenarios, @second_stage

# ---------------
# StochasticData
# ---------------

type StochasticData
    probability::Vector{Float64}
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
    # othervars::Vector{JuMP.Variable}
    othermap::Dict{JuMP.Variable,JuMP.Variable}
end

default_probability(m::JuMP.Model) = 1 / num_scenarios(m)
default_probability(::Void) = 1.0

# ---------------
# StochasticModel
# ---------------

# Constructor with the number of scenarios
function StochasticModel(;solver=JuMP.UnsetSolver(), parent=nothing, num_scenarios::Int=0, prob::Float64=default_probability(parent))
    m = JuMP.Model(solver=solver)
    if parent !== nothing
        stoch = getStochastic(parent)
        push!(stoch.children, m)
        push!(stoch.probability, prob)
    end
    m.ext[:Stochastic] = StochasticData(Float64[], JuMP.Model[], parent, num_scenarios, Dict{JuMP.Variable,JuMP.Variable}())
    m
end

# -------------
# Get functions
# -------------

getStochastic(m::JuMP.Model)  = m.ext[:Stochastic]
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
    return getProcIdxSet(numScens)
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
