# libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")

type UserData
    parent :: JuMP.Model
    child  :: JuMP.Model
end

cint(a::Int) = convert(Cint,a)
vcint(a::Vector{Int}) = convert(Vector{Cint},a)

function getConstraintTypes(m::JuMP.Model)
    numRows = length(m.linconstr)
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

function get_sparse_data(owner::JuMP.Model, interest::JuMP.Model, idx_set::Vector{Int})
    numRows = length(idx_set)
    rowptr   = Array(Int, numRows+1)

    # get a vague idea of how large submatrices will be
    nnz = 0
    for c in idx_set
        nnz += length(owner.linconstr[c].terms.coeffs)
    end

    colval   = Int[]
    rownzval = Float64[]

    nnz = 0
    tmprow   = JuMP.IndexedVector(Float64, owner.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for (it,c) in enumerate(idx_set)
        rowptr[it] = nnz + 1
        coeffs = owner.linconstr[c].terms.coeffs
        vars = owner.linconstr[c].terms.vars
        for (it,ind) in enumerate(coeffs)
            if vars[it].m == interest
                JuMP.addelt!(tmprow, vars[it].col, ind)
            end
        end
        for i in 1:tmprow.nnz
            nnz += 1
            idx = tmpnzidx[i]
            push!(colval, idx)
            push!(rownzval, tmpelts[idx])
        end
        JuMP.empty!(tmprow)
    end
    rowptr[numRows+1] = nnz + 1
    
    return rowptr, colval, rownzval
end

function get_sparse_Q(m::JuMP.Model)
    n = m.numCols
    vars1 = Int[x.col for x in m.obj.qvars1]
    vars2 = Int[x.col for x in m.obj.qvars2]
    coeff = m.obj.qcoeffs
    Q = sparse(vars1, vars2, coeff, n, n)
    istriu(Q) && (Q = Q')
    @assert istril(Q)
    return Q.colptr, Q.rowval, Q.nzval
end

include("pips_callbacks.jl")

function pips_solve(master::JuMP.Model)
    @assert getparent(master) == nothing # make sure this is master problem
    @assert length(getchildren(master)) == 1 # only one form of subproblem (at this point)

    child = getchildren(master)[1]

    # MPI data
    MPI.init()
    comm = MPI.COMM_WORLD

    root = 0
    size = MPI.size(comm)
    rank = MPI.rank(comm)

    numScens = num_scenarios(master)
    scenPerRank = iceil(numScens/size)

    user_data = UserData(master, child)

    eq_idx_m, ineq_idx_m = getConstraintTypes(master)
    eq_idx_c, ineq_idx_c = getConstraintTypes(child)

    n_eq_m, n_ineq_m = length(eq_idx_m), length(ineq_idx_m)
    n_eq_c, n_ineq_c = length(eq_idx_c), length(ineq_idx_c)

    val = ccall((libpips,"PIPSSolve"), Ptr{Void}, (Ptr{Cint},  # MPI_COMM
                                                   Ptr{Void},
                                                   Cint,       # numScens
                                                   Cint,       # nx0
                                                   Cint,       # my0
                                                   Cint,       # mz0
                                                   Cint,       # nx
                                                   Cint,       # my
                                                   Cint,       # mz
                                                   Any,  # Q
                                                   Any,  # nnzQ
                                                   Any,  # c
                                                   Any,  # A
                                                   Any,  # nnzA
                                                   Any,  # B
                                                   Any,  # nnzB
                                                   Any,  # b
                                                   Any,  # C
                                                   Any,  # nnzC
                                                   Any,  # D
                                                   Any,  # nnzD
                                                   Any,  # clow
                                                   Any,  # iclow
                                                   Any,  # cupp
                                                   Any,  # icupp
                                                   Any,  # xlow
                                                   Any,  # ixlow
                                                   Any,  # xupp
                                                   Any), # ixupp
                                                   &comm.fval,       
                                                   pointer_from_objref(user_data),
                                                   cint(numScens),   
                                                   cint(master.numCols),
                                                   cint(n_eq_m),        
                                                   cint(n_ineq_m),      
                                                   cint(child.numCols), 
                                                   cint(n_eq_c),        
                                                   cint(n_ineq_c),      
                                                   fQ,         
                                                   fnnzQ,      
                                                   fc,         
                                                   fA,         
                                                   fnnzA,      
                                                   fB,         
                                                   fnnZB,      
                                                   fb,         
                                                   fC,         
                                                   fnnzC,      
                                                   fD,         
                                                   fnnzD,      
                                                   fclow,      
                                                   ficlow,     
                                                   fcupp,      
                                                   ficupp,     
                                                   fxlow,      
                                                   fixlow,     
                                                   fxupp,      
                                                   fixupp)     
 
    MPI.finalize()
end
