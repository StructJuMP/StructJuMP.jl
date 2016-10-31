module StructJuMP

using JuMP # To reexport, should be using (not import)
import MathProgBase
import MathProgBase.MathProgSolverInterface
import ReverseDiffSparse

# These modules could be optional.
# import StructJuMPSolverInterface
# import MPI

export StructuredModel, getStructure, getparent, getchildren, getProcIdxSet,
       num_scenarios, @second_stage, getprobability, getMyRank
# Macro to exportall
macro exportall(pkg)
    Expr(:export, names(JuMP)...)
end
@exportall JuMP


# ---------------
# StructureData
# ---------------
if isdefined(:MPI)
    type MPIWrapper
        comm::MPI.Comm
        init::Function

        function MPIWrapper()
            instance = new(MPI.Comm(-1))
            finalizer(instance, freeMPIWrapper)

            instance.init = function(ucomm::MPI.Comm)
                if isdefined(:MPI) && MPI.Initialized() && ucomm.val == -1
                    instance.comm = MPI.COMM_WORLD
                elseif isdefined(:MPI) && !MPI.Initialized()
                    MPI.Init()
                    instance.comm = MPI.COMM_WORLD
                elseif isdefined(:MPI) && MPI.Initialized() && ucomm.val != -1
                    instance.comm = ucomm
                elseif isdefined(:MPI) && MPI.Finalized()
                    error("MPI is already finalized!")
                else
                    #doing nothing
                end
            end

            return instance
        end

    end
    function freeMPIWrapper(instance::MPIWrapper)
        if isdefined(:MPI) && MPI.Initialized() && !MPI.Finalized()
            MPI.Finalize()
        end
    end

    const mpiWrapper = MPIWrapper();

    type StructureData
        probability::Vector{Float64}
        children::Dict{Int,JuMP.Model}
        parent
        num_scen::Int
        othermap::Dict{JuMP.Variable,JuMP.Variable}
        mpiWrapper::MPIWrapper
    end
else
    type StructureData
        probability::Vector{Float64}
        children::Dict{Int,JuMP.Model}
        parent
        num_scen::Int
        othermap::Dict{JuMP.Variable,JuMP.Variable}
    end
end
default_probability(m::JuMP.Model) = 1 / num_scenarios(m)
default_probability(::Void) = 1.0

# ---------------
# StructuredModel
# ---------------

# Constructor with the number of scenarios
function StructuredModel(;solver=JuMP.UnsetSolver(), parent=nothing, same_children_as=nothing, id=0, comm=MPI.Comm(-1), num_scenarios::Int=0, prob::Float64=default_probability(parent))
    m = JuMP.Model(solver=solver)
    if parent === nothing
        id = 0
        if isdefined(:MPI)
            mpiWrapper.init(comm)
        end
        if isdefined(:StructJuMPSolverInterface)
            JuMP.setsolvehook(m,StructJuMPSolverInterface.sj_solve)
        end
    else
        @assert id != 0 
        stoch = getStructure(parent)
        stoch.children[id] = m
        push!(stoch.probability, prob)
    end

    if same_children_as !== nothing
      if !isa(same_children_as, JuMP.Model) || !haskey(same_children_as.ext, :Stochastic)
        error("The JuMP model given for the argument `same_children_as' is not valid. Please create it using the `StructuredModel' function.")
      end
      probability = same_children_as.ext[:Stochastic].probability
      children = same_children_as.ext[:Stochastic].children
    else
      probability = Float64[]
      children = Dict{Int, JuMP.Model}()
    end
    if isdefined(:MPI)
        m.ext[:Stochastic] = StructureData(probability, children, parent, num_scenarios, Dict{JuMP.Variable,JuMP.Variable}(), mpiWrapper)
    else
        m.ext[:Stochastic] = StructureData(probability, children, parent, num_scenarios, Dict{JuMP.Variable,JuMP.Variable}())
    end
    m
end

# -------------
# Get functions
# -------------

getStructure(m::JuMP.Model)   = m.ext[:Stochastic]::StructureData
getparent(m::JuMP.Model)      = getStructure(m).parent
getchildren(m::JuMP.Model)    = getStructure(m).children::Dict{Int,JuMP.Model}
getprobability(m::JuMP.Model) = getStructure(m).probability::Vector{Float64}
num_scenarios(m::JuMP.Model)  = getStructure(m).num_scen::Int


function getMyRank()
    myrank = 0;
    mysize = 1;
    if isdefined(:MPI) && MPI.Initialized() && !MPI.Finalized()
        comm = MPI.COMM_WORLD
        mysize = MPI.Comm_size(comm)
        myrank = MPI.Comm_rank(comm)
    end
    return myrank,mysize
end

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
