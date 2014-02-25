function passMasterSolution(m::Model)
    for ch in m.children
        ch.ext[:Stochastic].sol = m.colVal
    end
end

# Create sparse matrix A with all parent variables. Removes them from linear constraints afterwards.
# Once created, should be able to do A*x, where x is the solution to the master problem.
function createMasterMat(m::Model)
    stoch = getStochastic(m)
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


    stoch = getStochastic(m)
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
    stoch = getStochastic(m)
    @assert stoch.parent == nothing # make sure we're at the master problem

    master_stat solve(m)
    @assert master_stat == :Optimal
    passMasterSolution(m)

    for bl in stoch.children
        creatMasterMat(bl)
    end

    while true
        for bl in stoch.children
            solve(bl)
        end
        
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
