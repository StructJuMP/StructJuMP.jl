__precompile__()

module StructJuMP

using JuMP # To reexport, should be using (not import)
import MathProgBase
import ReverseDiffSparse

# These modules could be optional.
# import StructJuMPSolverInterface
# import MPI

export StructuredModel, getStructure, getparent, getchildren, getProcIdxSet,
       num_scenarios, @second_stage, getprobability, getMyRank,
       BendersBridge, DLP, loadAndSolveConicProblem
# Macro to exportall
macro exportall(pkg)
    Expr(:export, names(JuMP)...)
end
@exportall JuMP


# ---------------
# StructureData
# ---------------
type StructureData
    probability::Dict{Int,Float64}
    children::Dict{Int,JuMP.Model}
    parent
    num_scen::Int
    othermap::Dict{JuMP.Variable,JuMP.Variable}
    MPIWrapper # Empty unless StructJuMPwithMPI fills it
end
default_probability(m::JuMP.Model) = 1 / num_scenarios(m)
default_probability(::Void) = 1.0

# ---------------
# StructuredModel
# ---------------

function structprinthook(io::IO, m::Model)
    print(io, m, ignore_print_hook=true)
    print(io, "*** children ***\n")
    # TODO it would be nice to indent all the children
    # essentially wrap the IO object (subclass it) to add 4 spaces before each line
    # this would then recursively be wrapped as more stages are added
    for (id, mm) in getchildren(m)
      @printf(io, "Child ID %d:\n", id)
      print(io, mm)
      print(io, "\n")
    end
end

type DummyMPIWrapper
    comm::Int
    init::Function

    DummyMPIWrapper() = new(-1,identity)
end
const dummy_mpi_wrapper = DummyMPIWrapper()

# Constructor with the number of scenarios
function StructuredModel(;solver=JuMP.UnsetSolver(), parent=nothing, same_children_as=nothing, id=0, comm=nothing, num_scenarios::Int=0, prob::Float64=default_probability(parent), mpi_wrapper=dummy_mpi_wrapper)
    _comm = (comm == nothing ? mpi_wrapper.comm : comm)
    m = JuMP.Model(solver=solver)
    if parent === nothing
        id = 0
        mpi_wrapper.init(_comm)
        if isdefined(:StructJuMPSolverInterface)
            JuMP.setsolvehook(m,StructJuMPSolverInterface.sj_solve)
        end
    else
        @assert id != 0
        stoch = getStructure(parent)
        stoch.children[id] = m
        stoch.probability[id] = prob
    end

    if same_children_as !== nothing
        if !isa(same_children_as, JuMP.Model) || !haskey(same_children_as.ext, :Stochastic)
            error("The JuMP model given for the argument `same_children_as' is not valid. Please create it using the `StructuredModel' function.")
        end
        probability = same_children_as.ext[:Stochastic].probability
        children = same_children_as.ext[:Stochastic].children
    else
        probability = Dict{Int, Float64}()
        children = Dict{Int, JuMP.Model}()
    end
    m.ext[:Stochastic] = StructureData(probability, children, parent, num_scenarios, Dict{JuMP.Variable,JuMP.Variable}(), mpi_wrapper)

    # Printing children is important as well
    JuMP.setprinthook(m, structprinthook)

    m
end

# -------------
# Get functions
# -------------

getStructure(m::JuMP.Model)   = m.ext[:Stochastic]::StructureData
getparent(m::JuMP.Model)      = getStructure(m).parent
getchildren(m::JuMP.Model)    = getStructure(m).children::Dict{Int,JuMP.Model}
getprobability(m::JuMP.Model) = getStructure(m).probability::Dict{Int, Float64}
num_scenarios(m::JuMP.Model)  = getStructure(m).num_scen::Int

getProcIdxSet(dummy_mpi_wrapper::DummyMPIWrapper, num_scenarios) = 1:num_scenarios

function getProcIdxSet(m::JuMP.Model)
    haskey(m.ext[:Stochastic]) || error("Cannot use @second_stage without using the StructuredModel constructor")
    numScens = num_scenarios(m)
    return getProcIdxSet(getStructure(m).mpi_wrapper, numScens)
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

include("BendersBridge.jl")

end
