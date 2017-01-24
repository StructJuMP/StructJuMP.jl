function passMasterSolution(m::Model)
    for ch in m.children
        ch.ext[:Stochastic].sol = m.colVal
    end
end

# Create sparse matrix A with all parent variables. Removes them from linear constraints afterwards.
# Once created, should be able to do A*x, where x is the solution to the master problem.
function createMasterMat(m::Model)
    stoch = getStructure(m)
    parent = stoch.parent
    for (nrow,con) in enumerate(m.linconstr)
        aff = con.terms
        for (id,var) in zip(reverse(aff.vars), length(aff.vars):-1:1)
            if var.m == parent
                push!(I, nrow)
                push!(J, var.col)
                push!(V, aff.coeff[id])
                splice!(aff.vars, id)
                splice!(aff.coeffs, id)
            end
        end
    end
    stoch.parentMat = sparse(I,J,V,length(m.linconstr),parent.numCols)
end

function prepProblemBounds(m::Model)

    objaff::AffExpr = m.obj.aff
        
    # We already have dense column lower and upper bounds

    # Create dense objective vector
    f = zeros(m.numCols)
    for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] += objaff.coeffs[ind]
    end

    # Create row bounds
    numRows = length(m.linconstr)
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)


    stoch = getStructure(m)
    if stoch.parent != nothing && !isempty(stoch.parent.colVal)
        offset = stoch.parentMat * parent.colVal
        for c in 1:numRows
            rowlb[c] = m.linconstr[c].lb - offset[c]
            rowub[c] = m.linconstr[c].ub - offset[c]
        end
    else
        rowlb[c] = m.linconstr[c].lb
        rowub[c] = m.linconstr[c].ub
    end
    
    return f, rowlb, rowub
end

# Solve stochastic LP with Bender's decomposition
function solveStochastic(m::Model)
    stoch = getStructure(m)
    @assert stoch.parent == nothing # make sure we're at the master problem

    @variable(m, θ) # hopefully this name doesn't conflict...

    master_stat solve(m)
    @assert master_stat == :Optimal
    passMasterSolution(m)

    m.obj += θ

    for bl in stoch.children
        createMasterMat(bl)
    end

    while true
        for bl in stoch.children
            solve(bl)
            if bl.objVal > tol
            𝛑 = getconstrsolution(bl.internalModel)
            lhs = 𝛑 * getconstrmatrix(bl.internalModel)
            lb  = 𝛑 * getconstrLB(bl.internalModel) # what to do here...
            ub  = 𝛑 * getconstrUB(bl.internalModel) # and here...
            addconstr!(m.internalModel, Int64[1:length(lhs)], lhs, lb, ub)
        end
    end

end
