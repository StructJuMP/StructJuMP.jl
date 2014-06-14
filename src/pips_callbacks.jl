# TODO: cache sparsity information somewhere so we don't have to compute twice
function Q(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    rowptr, colvals, rownzvals = get_sparse_Q(host)
    unsafe_copy!(krowM, vcint(rowptr.-1), host.numCols+1)
    unsafe_copy!(jcolM, vcint(colvals.-1), length(colvals))
    unsafe_copy!(M, rownzvals, length(rownzvals))
    return nothing
end

function nnzQ(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, colvals, _ = get_sparse_Q(tar)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    return nothing
end

function A(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    eq_idx, _ = getConstraintTypes(host)
    rowptr, colvals, rownzvals = get_sparse_data(host, master, eq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     vcint(rownzvals), length(colvals))
    return nothing
end

function B(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    id == root && return nothing
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    eq_idx, _ = getConstraintTypes(child)
    rowptr, colvals, rownzvals = get_sparse_data(child, child, eq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     vcint(rownzvals), length(colvals))
    return nothing
end

function C(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, ineq_idx = getConstraintTypes(host)
    rowptr, colvals, rownzvals = get_sparse_data(host, master, ineq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     vcint(rownzvals), length(colvals))
    return nothing
end

function D(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    id == root && return nothing
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    _, ineq_idx = getConstraintTypes(child)
    rowptr, colvals, rownzvals = get_sparse_data(child, child, ineq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    n_eq+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     vcint(rownzvals), length(colvals))
    return nothing
end

function nnzA(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    eq_idx, _ = getConstraintTypes(host)
    _, colvals, _ = get_sparse_data(host, master, eq_idx)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    return nothing
end

function nnzB(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    if id == root
        unsafe_store!(nnz, cint(0), 1)
    else
        usr = unsafe_pointer_to_objref(user_data)::UserData
        master, child = usr.master, usr.child
        eq_idx, _ = getConstraintTypes(child)
        _, colvals, _ = get_sparse_data(child, child, eq_idx)
        unsafe_store!(nnz, cint(length(colvals)), 1)
    end
    return nothing
end

function nnzC(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, ineq_idx = getConstraintTypes(host)
    _, colvals, _ = get_sparse_data(host, master, ineq_idx)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    return nothing
end

function nnzD(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    if id == root
        unsafe_store!(nnz, cint(0), 1)
    else
        usr = unsafe_pointer_to_objref(user_data)::UserData
        master, child = usr.master, usr.child
        _, ineq_idx = getConstraintTypes(child)
        _, colvals, _ = get_sparse_data(child, child, ineq_idx)
        unsafe_store!(nnz, cint(length(colvals)), 1)
    end
    return nothing
end

function b(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    eq_idx, _ = getConstraintTypes(master)
    _, rlb, _ = JuMP.prepProblemBounds(master)
    @assert len == length(eq_idx)
    unsafe_copy!(vec, vcint(rlb[eq_idx]), len)
    return nothing
end

function c(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    f, _, _ = JuMP.prepProblemBounds(master)
    @assert len == length(host.numCols)
    unsafe_copy!(vec, vcint(f), len)
    return nothing
end

function clow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, rlb, _ = JuMP.prepProblemBounds(master)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(rlb[it]) ? 0.0 : rlb[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function cupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, _, rub = JuMP.prepProblemBounds(master)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(rlb[it]) ? 0.0 : rub[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function xlow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)

    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(host.colLower[it]) ? 0.0 : host.colLower[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function xupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(host.colUpper[it]) ? 0.0 : host.colUpper[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function iclow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, rlb, _ = JuMP.prepProblemBounds(master)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(rlb[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function icupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    _, _, rub = JuMP.prepProblemBounds(master)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(rlb[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function ixlow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)
    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(host.colLower[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
end

function ixupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master, child = usr.master, usr.child
    host = (id == root ? master : child)

    @assert len == length(rlb)
    for it in 1:len
        val = (isinf(host.colUpper[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
    end
    return nothing
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

for (name,src1,src2) in [(:iclow, :(get_ineq_idx(model)), :rlb_c),
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

for (mat_name,old_name) in [(:fQ,:Q),
                                     (:fA,:A),
                                     (:fB,:B),
                                     (:fC,:C),
                                     (:fD,:D)]
    @eval $mat_name =
            cfunction($old_name, Void, (Ptr{Void},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}))
end

for (nnz_name,old_name) in [(:fnnzQ,:Q),
                                     (:fnnzA,:A),
                                     (:fnnzB,:B),
                                     (:fnnzC,:C),
                                     (:fnnzD,:D)]
    @eval $nnz_name =
            cfunction($old_name, Void, (Ptr{Void},Cint,Ptr{Cint}))
end

for (vec_name,old_name) in [(:fb,:b), (:fc,:c),
                            (:fclow,:clow), (:fcupp,:cupp),
                            (:fxlow,:xlow), (:fxupp,:xupp),
                            (:ficlow,:iclow), (:ficupp,:icupp),
                            (:fixlow,:ixlow), (:fixupp,:ixupp)]
    @eval $vec_name =
        cfunction($old_name, Void, (Ptr{Void},Cint,Ptr{Cdouble},Cint))
end
