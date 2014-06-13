libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")

cint(a::Int) = convert(Cint,a)
vcint(a::Vector{Int}) = convert(Vector{Cint},a)

function get_sparse_data(owner::Model, interest::Model, idx_set::Vector{Int})
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
    for c in idx_set
        coeffs = owner.linconstr[c].terms.coeffs
        vars = owner.linconstr[c].terms.vars
        for (it,ind) in enumerate(coeffs)
            if vars[it].m == interest
                addelt!(tmprow, vars[it].col, coeffs[ind])
            end
        end
        for i in 1:tmprow.nnz
            nnz += 1
            idx = tmpnzidx[i]
            push!(colval, idx)
            push!(rownzval, tmpelts[idx])
        end
        empty!(tmprow)
    end
    rowptr[numRows+1] = nnz + 1
    
    return rowptr, colval, rownzval
end

function get_sparse_Q(m::Model)
    n = m.numCols
    qobj = m.obj
    vars1 = Int[x.col for x in m.obj.qvars1]
    vars2 = Int[x.col for x in m.obj.qvars2]
    coeff = m.obj.coeffs
    Q = sparse(vars1, vars2, coeff)
    @assert istriu(Q)
    return Q.rowptr, Q.colval, Q.nzval
end

function pips_solve(m::Model)
    @assert parent(m) == nothing # make sure this is master problem

    # MPI data
    comm = MPI.COMM_WORLD
    root = 0
    size = MPI.size(comm)
    rank = MPI.rank(comm)

    numScens = num_scenarios(m)
    scenPerRank = iceil(numScens/size)

    rank_to_model(id::Cint) = children(m)[rem1(id, scenPerRank)]

    f, rowlb, rowub = JuMP.prepProblemBounds(m)

    eq_idx, ineq_idx = getConstraintTypes(m)
    n_eq, n_ineq = length(eq_idx), length(ineq_idx)
    rlb, rub = rowlb[eq_idx], rowub[ineq_idx]

    #####################################################
    # Callback functions for matrices, vectors, and nnz's
    #####################################################

    # TODO: cache sparsity information somewhere so we don't have to compute twice
    function fQ(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        tar = rank_to_model(m,id)
        rowptr, colvals, rownzvals = get_sparse_Q(tar)
        unsafe_copy!(krowM, vcint(rowptr), tar.numCols+1)
        unsafe_copy!(jcolM, vcint(colvals), length(colvals))
        unsafe_copy!(M, rownzvals, length(rownzvals))
        return nothing
    end

    function fnnzQ(user_data, id::Cint, nnz::Ptr{Cint})
        tar = rank_to_model(m,id)
        rowptr, colvals, rownzvals = get_sparse_Q(tar)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        return nothing
    end

    for (name1, name2, typ) in [(:fA, :fnnzA, "eq"),
                                (:fB, :fnnzB, "eq"), 
                                (:fC, :fnnzC, "ineq"), 
                                (:fD, :fnnzD, "ineq")]
        @eval begin
            function $(name1)(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
                rowptr, colvals, rownzvals = get_sparse_data(m, rank_to_model(m,id), symbol(typ*"_idx"))
                unsafe_copy!(krowM, vcint(symbol(typ*"_rowptr[id]")),    n_eq+1)
                unsafe_copy!(jcolM, vcint(symbol(typ*"_colvals[id]")),   symbol("length("*typ*"_colvals[id])"))
                unsafe_copy!(M,     vcint(symbol(typ*"_rownzvals[id]")), symbol("length("*typ*"_colvals[id])"))
                return nothing
            end
            function $(name2)(user_data, id::Cint, nnz::Ptr{Cint})
                row_ptr, colvals, rownzvals = get_sparse_data(m, rank_to_model(m,id), symbol(typ*"_idx"))
                unsafe_store!(nnz, cint($(lngth)), 1)
                return nothing
            end
        end
    end

    for (name, lngth, stor) in [(:fb, :n_eq, :eq_rhs),
                                (:fc, :(length(f)), :f)]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                unsafe_copy!(vec, vcint($(stor)), len)
                return nothing
            end
        end
    end

    for (name1, name2, lngth, stor) in [(:fclow, :ficlow, :n_ineq, :rlb),
                                        (:fcupp, :ficupp, :n_ineq, :rub),
                                        (:fxlow, :fixlow, :(m.numCols), :(m.colLower)),
                                        (:fxupp, :fixupp, :(m.numCols), :(m.colUpper))]
        @eval begin
            function $(name1)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                for it in 1:len
                    val = (isinf($stor)[it]) ? 0.0 : $(stor)[it]
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
                return nothing
            end
            function $(name2)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                for it in 1:len
                    val = (isinf($(stor)[it]) ? 0.0 : 1.0)
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
                return nothing
            end
        end
    end

    ccall((libpips,:PIPSSolve), Void, (Ptr{Void},  # MPI_COMM
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
                                       Ptr{Void}), # ixupp
                                      (&comm.fval,       
                                       cint(numScens),   
                                       nx0,        
                                       my0,        
                                       mz0,        
                                       nx,         
                                       my,         
                                       mz,         
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
                                       fixupp))     

end
