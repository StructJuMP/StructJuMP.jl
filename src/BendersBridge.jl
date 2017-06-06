using JuMP
using StructJuMP

include("Benders_pmap.jl")

#===================================================
 The following function is taken from JuMP src/solvers.jl:conicdata(). 
 Modified for extension StructJuMP/Benders_pmap.jl.
===================================================#
function conicconstraintdata(m::Model)
    # TODO this generates redundant constraints for variables that are already bounded (one-sided)
    stoch = getStructure(m)
    parent = stoch.parent
    numMasterCols = 0
    if parent !== nothing
        numMasterCols = parent.numCols
    end
    v = Symbol[]

    var_cones = Any[cone for cone in m.varCones]
    con_cones = Any[]
    nnz = 0

    numSDPRows, numSymRows, nnz = JuMP.getSDrowsinfo(m)

    linconstr::Vector{LinearConstraint} = m.linconstr
    numLinRows = length(linconstr)

    numBounds = JuMP.variable_range_to_cone!(var_cones, m)

    nnz += numBounds
    for c in 1:numLinRows
        nnz += length(linconstr[c].terms.coeffs)
    end

    numSOCRows = JuMP.getNumSOCRows(m)
    numNormRows = length(m.socconstr)

    numRows = numLinRows + numBounds + numSOCRows + numSDPRows + numSymRows

    # constr_to_row is not used but fill_bounds_constr! and fillconstr! for SDP needs them
    constr_to_row = Array{Vector{Int}}(numBounds + 2*length(m.sdpconstr))

    b = Array{Float64}(numRows)

    I_m = Int[]
    J_m = Int[]
    V_m = Float64[]
    I_s = Int[]
    J_s = Int[]
    V_s = Float64[]

    # Fill it up
    if numMasterCols > 0
        tmprow_m = JuMP.IndexedVector(Float64, parent.numCols)
    end
    tmprow_s = JuMP.IndexedVector(Float64, m.numCols)

    JuMP.fillconstrRHS!(b, con_cones, 0, m.linconstr)
    if numMasterCols > 0
        JuMP.fillconstrLHS!(I_m, J_m, V_m, tmprow_m, 0, m.linconstr, parent, true)
    end
    c = JuMP.fillconstrLHS!(I_s, J_s, V_s, tmprow_s, 0, m.linconstr, m, true)

    for idx in 1:m.numCols
        # identify integrality information
        push!(v, m.colCat[idx])
    end
    c, d = JuMP.fill_bounds_constr!(I_s, J_s, V_s, b, con_cones, constr_to_row, c, 0, m)

    @assert c == numLinRows + numBounds
    @assert d == numBounds

    JuMP.fillconstrRHS!(b, con_cones, c, m.socconstr)
    if numMasterCols > 0
        JuMP.fillconstrLHS!(I_m, J_m, V_m, tmprow_m, c, m.socconstr, parent, true)
    end
    c = JuMP.fillconstrLHS!(I_s, J_s, V_s, tmprow_s, c, m.socconstr, m, true)

    @assert c == numLinRows + numBounds + numSOCRows

    if numMasterCols > 0
        c, d = JuMP.fillconstr!(I_m, J_m, V_m, b, con_cones, tmprow_m, constr_to_row, c, d, m.sdpconstr, m, true)
    end
    c, d = JuMP.fillconstr!(I_s, J_s, V_s, b, con_cones, tmprow_s, constr_to_row, c, d, m.sdpconstr, m, true)

    if c < length(b)
        # This happens for example when symmetry constraints are dropped with SDP
        resize!(b, c)
    end

    f_s = JuMP.prepAffObjective(m)

    # The conic MPB interface defines conic problems as
    # always being minimization problems, so flip if needed
    m.objSense == :Max && scale!(f_s, -1.0)

    if numMasterCols > 0
        JuMP.rescaleSDcols!(spzeros(numMasterCols), J_m, V_m, parent)
    end
    JuMP.rescaleSDcols!(f_s, J_s, V_s, m)

    A = sparse(I_m, J_m, V_m, numRows, numMasterCols)
    B = sparse(I_s, J_s, V_s, numRows, m.numCols)

    return f_s, A, B, b, var_cones, con_cones, v
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

    for i = 1:num_scenarios(m)
        (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(getchildren(m)[i])
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

    for i = 1:num_scenarios(m)
        (c,A,B,b,var_cones, constr_cones, v) = conicconstraintdata(getchildren(m)[i])
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
