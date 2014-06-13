libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")

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

    # rank_to_model(id::Cint) = getchildren(m)[rem1(id, scenPerRank)]

    f_m, rlb_m, rub_m = JuMP.prepProblemBounds(master)
    f_c, rlb_c, rub_c = JuMP.prepProblemBounds(child)

    eq_idx_m, ineq_idx_m = getConstraintTypes(master)
    eq_idx_c, ineq_idx_c = getConstraintTypes(child)

    n_eq_m, n_ineq_m = length(eq_idx_m), length(ineq_idx_m)
    n_eq_c, n_ineq_c = length(eq_idx_c), length(ineq_idx_c)

    eq_rhs_m = rlb_m[eq_idx_m]
    eq_rhs_c = rlb_c[eq_idx_c]

    ineq_lb_m, ineq_ub_m = rlb_m[ineq_idx_m], rub_m[ineq_idx_m]
    ineq_lb_c, ineq_ub_c = rlb_c[ineq_idx_c], rub_c[ineq_idx_c]

    #####################################################
    # Callback functions for matrices, vectors, and nnz's
    #####################################################

    # TODO: cache sparsity information somewhere so we don't have to compute twice
    function Q(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        host = (id == root ? master : child)
        rowptr, colvals, rownzvals = get_sparse_Q(host)
        unsafe_copy!(krowM, vcint(rowptr.-1), host.numCols+1)
        unsafe_copy!(jcolM, vcint(colvals.-1), length(colvals))
        unsafe_copy!(M, rownzvals, length(rownzvals))
        return nothing
    end

    function nnzQ(user_data, id::Cint, nnz::Ptr{Cint})
        host = (id == root ? master : child)
        _, colvals, _ = get_sparse_Q(tar)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        return nothing
    end

    function A(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        if id == root
            host = master
            rng = eq_idx_m
        else
            host = child
            rng = eq_idx_c
        end
        rowptr, colvals, rownzvals = get_sparse_data(host, master, rng)
        unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
        unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
        unsafe_copy!(M,     vcint(rownzvals), length(colvals))
        return nothing
    end

    function B(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        id == root && return nothing
        rowptr, colvals, rownzvals = get_sparse_data(child, child, eq_idx_c)
        unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
        unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
        unsafe_copy!(M,     vcint(rownzvals), length(colvals))
        return nothing
    end

    function C(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        if id == root
            host = master
            rng = ineq_idx_m
        else
            host = child
            rng = ineq_idx_c
        end
        rowptr, colvals, rownzvals = get_sparse_data(host, master, rng)
        unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
        unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
        unsafe_copy!(M,     vcint(rownzvals), length(colvals))
        return nothing
    end

    function D(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        id == root && return nothing
        rowptr, colvals, rownzvals = get_sparse_data(child, child, ineq_idx_m)
        unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
        unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
        unsafe_copy!(M,     vcint(rownzvals), length(colvals))
        return nothing
    end

    function nnzA(user_data, id::Cint, nnz::Ptr{Cint})
        if id == root
            host = master
            rng = eq_idx_m
        else
            host = child
            rng = eq_idx_c
        end
        _, colvals, _ = get_sparse_data(host, master, rng)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        return nothing
    end

    function nnzB(user_data, id::Cint, nnz::Ptr{Cint})
        if id == root
            unsafe_store!(nnz, cint(0), 1)
        else
            _, colvals, _ = get_sparse_data(child, child, eq_idx_c)
            unsafe_store!(nnz, cint(length(colvals)), 1)
        end
        return nothing
    end

    function nnzC(user_data, id::Cint, nnz::Ptr{Cint})
        if id == root
            host = master
            rng = eq_idx_m
        else
            host = child
            rng = eq_idx_c
        end
        _, colvals, _ = get_sparse_data(host, master, rng)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        return nothing
    end

    function nnzD(user_data, id::Cint, nnz::Ptr{Cint})
        if id == root
            unsafe_store!(nnz, cint(0), 1)
        else
            _, colvals, _ = get_sparse_data(child, child, ineq_idx_m)
            unsafe_store!(nnz, cint(length(colvals)), 1)
        end
        return nothing
    end

    function b(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
        if id == root
            @assert len == n_eq_m
            unsafe_copy!(vec, vcint(eq_rhs_m), len)
        else
            @assert len == n_eq_c
            unsafe_copy!(vec, vcint(eq_rhs_c), len)
        end
    end

    function c(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
        if id == root
            @assert len == n_eq_m
            unsafe_copy!(vec, vcint(f_m), len)
        else
            @assert len == n_eq_c
            unsafe_copy!(vec, vcint(f_c), len)
        end
    end

    for (name,src1,src2) in [(:clow, :rlb_m, :rlb_c),
                             (:cupp, :rub_m, :rub_c),
                             (:xlow, :(master.colLower), :(child.colUpper)),
                             (:xupp, :(master.colUpper), :(child.colUpper))]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                src = (id == root ? $src1 : $src2)
                @assert len == length(src)
                for it in 1:len
                    val = (isinf(src)[it]) ? 0.0 : src[it]
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
                return nothing
            end
        end
    end

    for (name,src1,src2) in [(:iclow, :rlb_m, :rlb_c),
                             (:icupp, :rub_m, :rub_c),
                             (:ixlow, :(master.colLower), :(child.colUpper)),
                             (:ixupp, :(master.colUpper), :(child.colUpper))]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                src = (id == root ? $src1 : $src2)
                @assert len == length(src)
                for it in 1:len
                    val = (isinf(src)[it]) ? 0.0 : 1.0
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
                return nothing
            end
        end
    end

    for name in [:Q, :A, :B, :C, :D]
        @eval begin 
            symbol("f"*name) = 
                cfunction($name, Void, (Ptr{Void},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}))
            symbol("fnnz"*name) = 
                cfunction($name, Void, (Ptr{Void},Cint,Cint,Ptr{Cint}))
        end
    end

    for name in [:b, :c, :clow, :cupp, :xlow, :xupp, :iclow, :icupp, :ixlow, :ixupp]
        @eval symbol("f"*name) = 
            cfunction($name, Void, (Ptr{Void},Cint,Ptr{Cdouble},Cint))
    end

    comm_val = Cint[comm.fval]

    val = ccall((libpips,"PIPSSolve"), Void, (Ptr{Cint},  # MPI_COMM
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
                                              comm_val,       
                                              cint(numScens),   
                                              cint(master.numCols),
                                              cint(n_eq_m),        
                                              cint(n_ineq_m),      
                                              cint(child.numCols), 
                                              cint(m_eq_c),        
                                              cint(m_ineq_c),      
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
