# This should not specify the file extension, because MacOSX has *.dylib.
# The library path should be specified in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH.
libpips = dlopen("libpipsipm-shared")
PIPSSolve = dlsym(libpips,:PIPSSolve)

type ProblemData
    Q :: SparseMatrixCSC{Cdouble,Cint}
    A :: SparseMatrixCSC{Cdouble,Cint}
    B :: SparseMatrixCSC{Cdouble,Cint}
    C :: SparseMatrixCSC{Cdouble,Cint}
    D :: SparseMatrixCSC{Cdouble,Cint}
    b :: Vector{Cdouble}
    c :: Vector{Cdouble}
    clow :: Vector{Cdouble}
    cupp :: Vector{Cdouble}
    xlow :: Vector{Cdouble}
    xupp :: Vector{Cdouble}
    iclow :: Vector{Cdouble}
    icupp :: Vector{Cdouble}
    ixlow :: Vector{Cdouble}
    ixupp :: Vector{Cdouble}
end

type UserData
    master   :: ProblemData
    children :: Dict{Cint,ProblemData}
end

import Base.convert
Base.convert(::Type{Ptr{Void}}, x::UserData) = x

const root = 0

cint(a::Int) = Base.convert(Cint,a)
vcint(a::Vector{Int}) = Base.convert(Vector{Cint},a)

function get_child_index(m::JuMP.Model, id)
    size = MPI.size(MPI.COMM_WORLD)
    scenPerRank = iceil(num_scenarios(m)/size)
    rem1(id+1,scenPerRank)
end

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
    sizehint(colval, nnz)
    rownzval = Float64[]
    sizehint(rownzval, nnz)

    nnz = 0
    tmprow   = JuMP.IndexedVector(Float64, interest.numCols)
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

    mat = SparseMatrixCSC(interest.numCols, numRows, vcint(rowptr), vcint(colval), rownzval)
    mat = (mat')' # ugly ugly ugly
    mat.colptr .-= 1
    mat.rowval .-= 1
    return mat
end

function get_sparse_Q(m::JuMP.Model)
    n = m.numCols
    vars1 = Int[x.col for x in m.obj.qvars1]
    vars2 = Int[x.col for x in m.obj.qvars2]
    coeff_copy = copy(m.obj.qcoeffs)
    for i in 1:length(vars1)
        if vars1[i] == vars2[i] # "terms" form
            coeff_copy[i] *= 2
        end
        if vars1[i] > vars2[i]
            vars1[i], vars2[i] = vars2[i], vars1[i]
        end
    end
    Q = sparse(vcint(vars1), vcint(vars2), coeff_copy, n, n)
    istriu(Q) && (Q = Q')
    @assert istril(Q)
    Q.colptr .-= 1 # zero-indexed counting
    Q.rowval .-= 1
    return Q
end

include("pips_callbacks.jl")

function pips_solve(master::JuMP.Model)
    tic()
    @assert getparent(master) == nothing # make sure this is master problem
    @assert master.objSense == :Min

    children = getchildren(master)

    passToPIPS = !(length(children) == 0)
    child    = passToPIPS ? children[1] : master

    # MPI data
    comm = MPI.COMM_WORLD

    size = MPI.Comm_size(comm)
    rank = MPI.Comm_rank(comm)

    numScens = num_scenarios(master)
    scenPerRank = iceil(numScens/size)

    eq_idx_m, ineq_idx_m = getConstraintTypes(master)
    eq_idx_c, ineq_idx_c = getConstraintTypes(child) # should probably check this is const across children

    n_eq_m, n_ineq_m = length(eq_idx_m), length(ineq_idx_m)
    n_eq_c, n_ineq_c = length(eq_idx_c), length(ineq_idx_c)

    f, rlb, rub = JuMP.prepProblemBounds(master)
    Q = get_sparse_Q(master)
    A = get_sparse_data(master, master, eq_idx_m)
    B = SparseMatrixCSC(n_eq_m, master.numCols, fill(cint(1),n_eq_m+1),Cint[],Cdouble[])
    B.colptr .-= 1
    C = get_sparse_data(master, master, ineq_idx_m)
    D = SparseMatrixCSC(n_ineq_m, master.numCols, fill(cint(1),n_ineq_m+1),Cint[],Cdouble[])
    D.colptr .-= 1
    b = rlb[eq_idx_m]
    c = f
    clow = rlb[ineq_idx_m]
    cupp = rub[ineq_idx_m]
    iclow = Array(Cdouble, n_ineq_m)
    icupp = Array(Cdouble, n_ineq_m)
    for (it,val) in enumerate(ineq_idx_m)
        clow[it] = (isinf(rlb[val]) ? 0.0 : rlb[val])
        cupp[it] = (isinf(rub[val]) ? 0.0 : rub[val])
        iclow[it] = (isinf(rlb[val]) ? 0.0 : 1.0)
        icupp[it] = (isinf(rub[val]) ? 0.0 : 1.0)
    end
    xlow = Array(Cdouble, master.numCols)
    xupp = Array(Cdouble, master.numCols)
    ixlow = Array(Cdouble, master.numCols)
    ixupp = Array(Cdouble, master.numCols)
    for it in 1:master.numCols
        xlow[it] = (isinf(master.colLower[it]) ? 0.0 : master.colLower[it])
        xupp[it] = (isinf(master.colUpper[it]) ? 0.0 : master.colUpper[it])
        ixlow[it] = (isinf(master.colLower[it]) ? 0.0 : 1.0)
        ixupp[it] = (isinf(master.colUpper[it]) ? 0.0 : 1.0)
    end
    master_prob = ProblemData(Q,A,B,C,D,b,c,clow,cupp,xlow,xupp,iclow,icupp,ixlow,ixupp)

    #children_probs = Array(ProblemData, length(children))
    children_probs = Dict{Cint,ProblemData}()
    for (idx,child) in enumerate(children)
        f, rlb, rub = JuMP.prepProblemBounds(child)
        Q = get_sparse_Q(child)
        A = get_sparse_data(child, master, eq_idx_c)
        B = get_sparse_data(child, child, eq_idx_c)
        C = get_sparse_data(child, master, ineq_idx_c)
        D = get_sparse_data(child, child, ineq_idx_c)
        b = rlb[eq_idx_c]
        c = f
        clow = rlb[ineq_idx_c]
        cupp = rub[ineq_idx_c]
        iclow = Array(Cdouble, n_ineq_c)
        icupp = Array(Cdouble, n_ineq_c)
        for (it,val) in enumerate(ineq_idx_c)
            clow[it] = (isinf(rlb[val]) ? 0.0 : rlb[val])
            cupp[it] = (isinf(rub[val]) ? 0.0 : rub[val])
            iclow[it] = (isinf(rlb[val]) ? 0.0 : 1.0)
            icupp[it] = (isinf(rub[val]) ? 0.0 : 1.0)
        end
        xlow = Array(Cdouble, child.numCols)
        xupp = Array(Cdouble, child.numCols)
        ixlow = Array(Cdouble, child.numCols)
        ixupp = Array(Cdouble, child.numCols)
        for it in 1:child.numCols
            xlow[it] = (isinf(child.colLower[it]) ? 0.0 : child.colLower[it])
            xupp[it] = (isinf(child.colUpper[it]) ? 0.0 : child.colUpper[it])
            ixlow[it] = (isinf(child.colLower[it]) ? 0.0 : 1.0)
            ixupp[it] = (isinf(child.colUpper[it]) ? 0.0 : 1.0)
        end

        ind = convert(Cint, scenPerRank*rank + idx)
        children_probs[ind] = ProblemData(Q,A,B,C,D,b,c,clow,cupp,xlow,xupp,iclow,icupp,ixlow,ixupp)
    end

    gc_disable() # probably very dangerous...
    user_data = UserData(master_prob, children_probs)

    obj_val = [0.0]
    first_primal  = Array(Cdouble, master.numCols)
    second_primal = Array(Cdouble, numScens*child.numCols)
    first_dual    = Array(Cdouble, n_eq_m+n_ineq_m)
    second_dual   = Array(Cdouble, numScens*(n_eq_c+n_ineq_c))

    MPI.Barrier(comm)
    t1 = toc()
    rank == 0 && println("jump time #2: $t1 secs")

    if !passToPIPS
        # TODO: call Q, nnzQ, etc. to precompile
        gc_enable()
        return t1, 0.0
    end

    tic()
    ccall(PIPSSolve, Void, (Ptr{Cint},  # MPI_COMM
    #val = ccall(("PIPSSolve",libpips), Void, (Ptr{Void},  # MPI_COMM
                                                   Ptr{Void},
                                                   Cint,       # numScens
                                                   Cint,       # nx0
                                                   Cint,       # my0
                                                   Cint,       # mz0
                                                   Cint,       # nx
                                                   Cint,       # my
                                                   Cint,       # mz
                                                   Ptr{Void},  # Q
                                                   Ptr{Void},  # nnzQ
                                                   Ptr{Void},  # c
                                                   Ptr{Void},  # A
                                                   Ptr{Void},  # nnzA
                                                   Ptr{Void},  # B
                                                   Ptr{Void},  # nnzB
                                                   Ptr{Void},  # b
                                                   Ptr{Void},  # C
                                                   Ptr{Void},  # nnzC
                                                   Ptr{Void},  # D
                                                   Ptr{Void},  # nnzD
                                                   Ptr{Void},  # clow
                                                   Ptr{Void},  # iclow
                                                   Ptr{Void},  # cupp
                                                   Ptr{Void},  # icupp
                                                   Ptr{Void},  # xlow
                                                   Ptr{Void},  # ixlow
                                                   Ptr{Void},  # xupp
                                                   Ptr{Void}, # ixupp
                                                   Ptr{Cdouble}, # objval
                                                   Ptr{Cdouble}, # first-stage primal sol
                                                   Ptr{Cdouble}, # second-stage primal so
                                                   Ptr{Cdouble}, # first-stage dual
                                                   Ptr{Cdouble}), #second-stage dual
                                                   Cint[root],
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
                                                   fnnzB,
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
                                                   fixupp,
                                                   obj_val,
                                                   first_primal,
                                                   second_primal,
                                                   first_dual,
                                                   second_dual)
    gc_enable()
    MPI.Barrier(comm)
    t2 = toc()
    #MPI.finalize()
    return t1, t2
end
