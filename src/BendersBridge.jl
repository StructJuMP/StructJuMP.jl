using JuMP
using StructJuMP

include("Benders_pmap.jl")

#===================================================
 The following function is taken from JuMP
 src/solvers.jl. Modified for extension 
 StructJuMP/Benders_pmap.jl.
===================================================#
function conicconstraintdata(m::Model)
    # TODO this generates redundant constraints for variables that are already bounded (one-sided)
    stoch = getStructure(m)
    parent = stoch.parent
    numMasterCols = 0
    if parent != nothing
        numMasterCols = parent.numCols
    end
    v = Symbol[]
    obj_coeff = zeros(m.numCols)
    for i in eachindex(m.obj.aff.vars)
        var = m.obj.aff.vars[i]
        coeff = m.obj.aff.coeffs[i]
        obj_coeff[var.col] = coeff
    end

    var_cones = Any[cone for cone in m.varCones]
    con_cones = Any[]
    nnz = 0

    linconstr::Vector{LinearConstraint} = m.linconstr
    numLinRows = length(linconstr)
    numBounds = 0
    nonNeg  = Int[]
    nonPos  = Int[]
    free    = Int[]
    zeroVar = Int[]
    for i in 1:m.numCols
        seen = false
        lb, ub = m.colLower[i], m.colUpper[i]
        for (_,cone) in m.varCones
            if i in cone
                seen = true
                @assert lb == -Inf && ub == Inf
                break
            end
        end

        if !seen
            if lb != -Inf
                numBounds += 1
            end
            if ub != Inf
                numBounds += 1
            end
            if lb == 0 && ub == 0
                push!(zeroVar, i)
            elseif lb == 0
                push!(nonNeg, i)
            elseif ub == 0
                push!(nonPos, i)
            else
                push!(free, i)
            end
        end
    end

    if !isempty(zeroVar)
        push!(var_cones, (:Zero,zeroVar))
    end
    if !isempty(nonNeg)
        push!(var_cones, (:NonNeg,nonNeg))
    end
    if !isempty(nonPos)
        push!(var_cones, (:NonPos,nonPos))
    end
    if !isempty(free)
        push!(var_cones, (:Free,free))
    end

    nnz += numBounds
    for c in 1:numLinRows
        nnz += length(linconstr[c].terms.coeffs)
    end

    numSOCRows = 0
    for con in m.socconstr
        numSOCRows += length(con.normexpr.norm.terms) + 1
    end
    numRows = numLinRows + numBounds + numSOCRows

    b = Array(Float64, numRows)

    I_m = Int[]
    J_m = Int[]
    V_m = Float64[]
    I_s = Int[]
    J_s = Int[]
    V_s = Float64[]

    # Fill it up
    nnz = 0
    tmprow = JuMP.IndexedVector(Float64,m.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    nonneg_rows = Int[]
    nonpos_rows = Int[]
    eq_rows     = Int[]
    for c in 1:numLinRows
        if linconstr[c].lb == -Inf
            b[c] = linconstr[c].ub
            push!(nonneg_rows, c)
        elseif linconstr[c].ub == Inf
            b[c] = linconstr[c].lb
            push!(nonpos_rows, c)
        elseif linconstr[c].lb == linconstr[c].ub
            b[c] = linconstr[c].lb
            push!(eq_rows, c)
        else
            error("We currently do not support ranged constraints with conic solvers")
        end

        JuMP.assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars = linconstr[c].terms.vars
        # eliminated collect duplicates
        for ind in eachindex(coeffs)
            if vars[ind].m === parent
                push!(I_m, c)
                push!(J_m, vars[ind].col)
                push!(V_m, coeffs[ind])
            else
                push!(I_s, c)
                push!(J_s, vars[ind].col)
                push!(V_s, coeffs[ind])
            end
        end
    end

    c = numLinRows
    bndidx = 0
    for idx in 1:m.numCols
        lb = m.colLower[idx]
        # identify integrality information
        push!(v, m.colCat[idx])
        if lb != -Inf
            bndidx += 1
            nnz += 1
            c   += 1
            push!(I_s, c)
            push!(J_s, idx)
            push!(V_s, 1.0)
            b[c] = lb
            push!(nonpos_rows, c)
        end
        ub = m.colUpper[idx]
        if ub != Inf
            bndidx += 1
            c   += 1
            push!(I_s, c)
            push!(J_s, idx)
            push!(V_s, 1.0)
            b[c] = ub
            push!(nonneg_rows, c)
        end
    end

    if !isempty(nonneg_rows)
        push!(con_cones, (:NonNeg,nonneg_rows))
    end
    if !isempty(nonpos_rows)
        push!(con_cones, (:NonPos,nonpos_rows))
    end
    if !isempty(eq_rows)
        push!(con_cones, (:Zero,eq_rows))
    end
    @assert c == numLinRows + numBounds

    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    socidx = 0
    for con in m.socconstr
        socidx += 1
        expr = con.normexpr
        c += 1
        soc_start = c
        JuMP.collect_expr!(m, tmprow, expr.aff)
        nnz = tmprow.nnz
        indices = tmpnzidx[1:nnz]
        vars = expr.aff.vars
        for i in eachindex(vars)
            if vars[i].m === parent
                push!(I_m, c)
                push!(J_m, indices[i])
                push!(V_m, tmpelts[indices[i]])
            else
                push!(I_s, c)
                push!(J_s, indices[i])
                push!(V_s, tmpelts[indices[i]])
            end
        end
        b[c] = -expr.aff.constant
        for term in expr.norm.terms
            c += 1
            JuMP.collect_expr!(m, tmprow, term)
            nnz = tmprow.nnz
            indices = tmpnzidx[1:nnz]
            vars = term.vars
            for i = 1:length(vars)
                if vars[i].m == parent
                    push!(I_m, c)
                    push!(J_m, indices[i])
                    push!(V_m, -expr.coeff*tmpelts[indices[i]])
                else
                    push!(I_s, c)
                    push!(J_s, indices[i])
                    push!(V_s, -expr.coeff*tmpelts[indices[i]])
                end
            end
            b[c] = expr.coeff*term.constant
        end
        push!(con_cones, (:SOC, soc_start:c))
    end
    @assert c == numLinRows + numBounds + numSOCRows

    A = sparse(I_m, J_m, V_m, numRows, numMasterCols)
    B = sparse(I_s, J_s, V_s, numRows, m.numCols)
  
    return obj_coeff, A, B, b, var_cones, con_cones, v
end

function BendersBridge(m::Model, master_solver, sub_solver)

    c_all = Vector{Float64}[]
    A_all = SparseMatrixCSC{Float64}[]
    B_all = SparseMatrixCSC{Float64}[]
    b_all = Vector{Float64}[]
    K_all = Any[]
    C_all = Any[]
    v_all = Symbol[]

    (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(m)
    push!(c_all, c)
    push!(A_all, B)
    push!(b_all, b)
    push!(K_all, constr_cones)
    push!(C_all, var_cones)
    append!(v_all, v)

    for i = 1:length(m.ext[:Stochastic].children)
        (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(m.ext[:Stochastic].children[i])
        push!(c_all, getprobability(m)[i] * c)
        push!(A_all, A)
        push!(B_all, B)
        push!(b_all, b)
        push!(K_all, constr_cones)
        push!(C_all, var_cones)
        append!(v_all, v)
    end

    return Benders_pmap(c_all,A_all,B_all,b_all,K_all,C_all,v_all,master_solver,sub_solver)
end

function DLP(m::Model, solver)

    c_all = Vector{Float64}[]
    A_all = SparseMatrixCSC{Float64}[]
    B_all = SparseMatrixCSC{Float64}[]
    b_all = Vector{Float64}[]
    K_all = Any[]
    C_all = Any[]
    v_all = Symbol[]

    (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(m)
    push!(c_all, c)
    push!(A_all, B)
    push!(b_all, b)
    push!(K_all, constr_cones)
    push!(C_all, var_cones)
    append!(v_all, v)

    for i = 1:length(m.ext[:Stochastic].children)
        (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(m.ext[:Stochastic].children[i])
        push!(c_all, getprobability(m)[i] * c)
        push!(A_all, A)
        push!(B_all, B)
        push!(b_all, b)
        push!(K_all, constr_cones)
        push!(C_all, var_cones)
        append!(v_all, v)
    end
    # v_all is unused for now

    dlp_rows = sum(length(b) for b in b_all)
    dlp_cols = sum(length(c) for c in c_all)

    c_dlp = zeros(dlp_cols)
    A_dlp = spzeros(dlp_rows, dlp_cols)
    b_dlp = zeros(dlp_rows)
    K_dlp = Any[] # constr_cones
    C_dlp = Any[] # var_cones

    rows = 1
    cols = 1
    for i = 1:length(c_all)
        rows_end = rows + length(b_all[i]) - 1
        cols_end = cols + length(c_all[i]) - 1

        c_dlp[cols:cols_end] = c_all[i]
        b_dlp[rows:rows_end] = b_all[i]
        
        A_dlp[rows:rows_end, 1:size(A_all[i], 2)] = A_all[i]
        if i != 1
            A_dlp[rows:rows_end, cols:cols_end] = B_all[i-1]
        end
        
        append!(K_dlp, [(cone, rows-1 + idx) for (cone, idx) in K_all[i]])
        append!(C_dlp, [(cone, cols-1 + idx) for (cone, idx) in C_all[i]])
        rows = rows_end+1
        cols = cols_end+1
    end
    
    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c_dlp, A_dlp, b_dlp, K_dlp, C_dlp)
 
    # solve conic model
    MathProgBase.optimize!(model)
    status = MathProgBase.status(model)

    # return status, objective value, and solution
    return status, MathProgBase.getobjval(model), MathProgBase.getsolution(model)[1:length(c_all[1])]

end
