module JuMPStoch

# import JuMP.JuMPDict
importall JuMP

using MathProgBase
using MathProgBase.MathProgSolverInterface

importall Base

export StochasticData, StochasticModel, getStochastic, StochasticBlock

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

pushchild!(m::Model, block) = push!(m.ext[:Stochastic].children, block)

type StochasticData
    id
    children::Vector{Model}
    parent
end

StochasticData() = StochasticData(nothing,Model[],nothing)

function StochasticModel(;solver=nothing)
    m = Model(solver=solver)
    m.ext[:Stochastic] = StochasticData()
    return m
end

function StochasticModel(id, children, parent)
    m = Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(id, children, parent)
    return m
end

function getStochastic(m::Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

function StochasticBlock(m::Model, id)
    stoch = getStochastic(m)
    ch = StochasticModel(id, Model[], m)
    pushchild!(m, ch)
    return ch
end

function solveStochastic(m::Model)
    f, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)
    m.internalModel = model(m.solver)
    loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)

    for bl in m.blocks
        f, rowlb, rowub = prepProblemBounds(bl)
        A = prepConstrMatrix(bl)
        loadproblem!(bl.internalModel, A, bl.colLower, bl.colUpper, f, rowlb, rowub, bl.objSense)
    end

    # ape JuliaBenders
    nscen = length(m.blocks[1].children)
    converged = false
    niter = 0
    mastertime = 0.

    while true
        Tx = d.Tmat*stage1sol
        # solve benders subproblems
        nviolated = 0
        for s in 1:nscen
            optval, subgrad = solveSubproblem(scenarioData[s][1]-Tx,scenarioData[s][2]-Tx)
            if (optval > thetasol[s] + 1e-7)
                nviolated += 1
                addCut(clpmaster, optval, subgrad, stage1sol, s)
            end

        end

        if nviolated == 0
            break
        end
        println("Generated $nviolated violated cuts")
        # resolve master
        t = time()
        initial_solve(clpmaster)
        mastertime += time() - t
        @assert is_proven_optimal(clpmaster)
        sol = get_col_solution(clpmaster)
        stage1sol = sol[1:ncol1]
        thetasol = sol[(ncol1+1):end]
        niter += 1
    end

end

end
