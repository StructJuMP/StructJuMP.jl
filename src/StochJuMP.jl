module StochJuMP

import JuMP.JuMPDict
import JuMP.@gendict
using JuMP

import MPI

using MathProgBase
using MathProgBase.MathProgSolverInterface

import Base.parent
importall Base
using Base.Meta

export StochasticData, StochasticModel, getStochastic, parent, children, StochasticBlock

# JuMP rexports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint,
# Functions
    # Relevant to all
    print,show,
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addVar, addVars, solve, copy,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    # Macros and support functions
    @addConstraint, @defVar,
    @defConstrRef, @setObjective, addToExpression

type StochasticData
    id
    children::Vector{Model}
    parent
    num_scen::Int
end

StochasticData() = StochasticData(nothing,Model[],nothing,0)

function StochasticModel(;solver=nothing)
    m = Model(solver=solver)
    m.ext[:Stochastic] = StochasticData()
    return m
end

function StochasticModel(id, children, parent)
    m = Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(id, children, parent,0)
    return m
end

function getStochastic(m::Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

parent(m::Model) = getStochastic(m).parent
children(m::Model) = getStochastic(m).children
num_scenarios(m::Model) = getStochastic(m).num_scen

function StochasticBlock(m::Model, id)
    stoch = getStochastic(m)
    ch = StochasticModel(id, Model[], m)
    push!(getStochastic(m).children, ch)
    return ch
end

# function fill_sparse_data(m::Model, idx_set::Vector{Int})
#     numRows = length(idx_set)

#     # get a vague idea of how large submatrices will be
#     nnz = 0
#     for c in idx_set
#         nnz += length(m.linconstr[c].terms.coeffs)
#     end

#     nnzs      = Dict{JuMP.Model, Int}()
#     rowptrs   = Dict{JuMP.Model, Vector{Int}}()
#     colvals   = Dict{JuMP.Model, Vector{Int}}()
#     rownzvals = Dict{JuMP.Model, Vector{Float64}}()
#     tmprows   = Dict{JuMP.Model, JuMP.IndexedVector}()
#     for anc in ancestors
#         nnzs[anc] = 0
#         rowptrs[anc]   = Array(Int,num+1)
#         colvals[anc] = Int[]
#         sizehint(colvals[anc], nnz)
#         eq_rownzvals[anc] = Float64[]
#         sizehint(rownzvals[anc], nnz)
#         tmprows[anc] = IndexedVector(Float64, anc.numCols)
#     end

#     for c in idx_set
#         coeffs = m.linconstr[c].terms.coeffs
#         vars = m.linconstr[c].terms.vars
#         for (it,ind) in enumerate(coeffs)
#             addelt!(tmprows[vars[it].m], vars[it].m, coeffs[ind])
#         end
#         for anc in ancestors
#             tmprow = tmprows[anc]
#             for i in 1:tmprow.nnz
#                 nnzs[anc] += 1
#                 idx = tmprow.nzidx[i]
#                 push!(colvals[anc], idx)
#                 push!(rownzvals[anc], tmprow.elts[idx])
#             end
#         end
#         map(empty!, tmprows)
#     end

#     mats = Array(SparseMatrixCSC, length(ancestors))
#     for (it,anc) in enumerate(ancestors)
#         rowptrs[anc][num+1]   = nnzs[anc] + 1
#         mats[it] = SparseMatrixCSC(anc.numCols,
#                                     numRows,
#                                     rowptrs[anc],
#                                     colval[anc],
#                                     rownzval[anc])
#     end
#     return mats
# end

function getConstraintTypes(m::Model)
    eq_idx   = Int[]
    sizehint(eq_idx, numRows)
    ineq_idx = Int[]
    sizehint(ineq_idx, numRows)
    for it in 1:numRows
        if m.linconstr[it].lb == m.linconstr[it].ub
            push!(eq_idx, it)
        else
            push!(ineq_idx, it)
        end
    end
    return eq_idx, ineq_idx
end

# # ancestors[1] = current model
# constructMatrices(m::Model) = constructMatrices(Model[m])
# function constructMatrices(ancestors::Vector{Model})
#     m = ancestors[1]
#     # determine number of inequalities and equalities
#     eq_idx, ineq_idx = getConstraintTypes(m)

#     eq_mats   = fill_sparse_data(m, eq_idx)
#     ineq_mats = fill_sparse_data(m, ineq_idx)

#     eq_rhs  = m.colLower[eq_idx]
#     ineq_lb = m.colLower[ineq_idx]
#     ineq_ub = m.colUpper[ineq_idx]

#     return eq_mats, eq_rhs, ineq_mats, ineq_lb, ineq_ub
# end

# function prepStochasticConstrMatrix(m::Model)

#     stoch = getStochastic(m)
#     A0, b0, C0, d0_l, d0_u  = constructMatrices(m)

#     for child in children(m)
#         A_i, b_i, C_i, di_l, di_u = constructMatrices([m, child])
#     end

# end

include("pips.jl")

end
