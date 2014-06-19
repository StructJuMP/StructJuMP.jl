# TODO: cache sparsity information somewhere so we don't have to compute twice

function Q(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    rowptr, colvals, rownzvals = get_sparse_Q(host)
    unsafe_copy!(krowM, vcint(rowptr.-1), host.numCols+1)
    unsafe_copy!(jcolM, vcint(colvals.-1), length(colvals))
    unsafe_copy!(M, pointer(rownzvals), length(rownzvals))
    println("Q ($id):")
    println("rowptr    = $rowptr")
    println("colvals   = $colvals")
    println("rownzvals = $rownzvals")
    println()
    return nothing
end

function nnzQ(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    _, colvals, _ = get_sparse_Q(host)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    return nothing
end

function A(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master = usr.master
    child_id = get_child_index(master, id)
    host = (id == root ? master : usr.children[child_id])
    eq_idx, _ = getConstraintTypes(host)
    rowptr, colvals, rownzvals = get_sparse_data(host, master, eq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    length(eq_idx)+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     pointer(rownzvals), length(colvals))
    println("A ($id):")
    println("rowptr    = $rowptr")
    println("colvals   = $colvals")
    println("rownzvals = $rownzvals")
    println()
    return nothing
end

function B(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    id == root && return nothing
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    child = usr.children[child_id]
    eq_idx, _ = getConstraintTypes(child)
    rowptr, colvals, rownzvals = get_sparse_data(child, child, eq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    length(eq_idx)+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     pointer(rownzvals), length(colvals))
    println("B ($id):")
    println("rowptr    = $rowptr")
    println("colvals   = $colvals")
    println("rownzvals = $rownzvals")
    println()
    return nothing
end

function C(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master = usr.master
    child_id = get_child_index(master, id)
    host = (id == root ? master : usr.children[child_id])
    _, ineq_idx = getConstraintTypes(host)
    rowptr, colvals, rownzvals = get_sparse_data(host, master, ineq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    length(ineq_idx)+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     pointer(rownzvals), length(colvals))
    println("C ($id):")
    println("rowptr    = $rowptr")
    println("colvals   = $colvals")
    println("rownzvals = $rownzvals")
    println()
    return nothing
end

function D(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
    id == root && return nothing
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    child = usr.children[child_id]
    _, ineq_idx = getConstraintTypes(child)
    rowptr, colvals, rownzvals = get_sparse_data(child, child, ineq_idx)
    unsafe_copy!(krowM, vcint(rowptr.-1),    length(ineq_idx)+1)
    unsafe_copy!(jcolM, vcint(colvals.-1),   length(colvals))
    unsafe_copy!(M,     pointer(rownzvals), length(colvals))
    println("D ($id):")
    println("rowptr    = $rowptr")
    println("colvals   = $colvals")
    println("rownzvals = $rownzvals")
    println()
    return nothing
end

function nnzA(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master = usr.master
    child_id = get_child_index(master, id)
    host = (id == root ? master : usr.children[child_id])
    eq_idx, _ = getConstraintTypes(host)
    _, colvals, _ = get_sparse_data(host, master, eq_idx)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    println("nnzA ($id) = $(length(colvals))")
    return nothing
end

function nnzB(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    if id == root
        unsafe_store!(nnz, cint(0), 1)
        println("nnzB ($id) = 0")
    else
        usr = unsafe_pointer_to_objref(user_data)::UserData
        child_id = get_child_index(usr.master, id)
        child = usr.children[child_id]
        eq_idx, _ = getConstraintTypes(child)
        _, colvals, _ = get_sparse_data(child, child, eq_idx)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        println("nnzB ($id) = $(length(colvals))")
    end
    return nothing
end

function nnzC(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    usr = unsafe_pointer_to_objref(user_data)::UserData
    master = usr.master
    child_id = get_child_index(master, id)
    host = (id == root ? master : usr.children[child_id])
    _, ineq_idx = getConstraintTypes(host)
    _, colvals, _ = get_sparse_data(host, master, ineq_idx)
    unsafe_store!(nnz, cint(length(colvals)), 1)
    println("nnzC ($id) = $(length(colvals))")
    return nothing
end

function nnzD(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
    if id == root
        unsafe_store!(nnz, cint(0), 1)
        println("nnzD ($id) = 0")
    else
        usr = unsafe_pointer_to_objref(user_data)::UserData
        child_id = get_child_index(usr.master, id)
        child = usr.children[child_id]
        _, ineq_idx = getConstraintTypes(child)
        _, colvals, _ = get_sparse_data(child, child, ineq_idx)
        unsafe_store!(nnz, cint(length(colvals)), 1)
        println("nnzD ($id) = $(length(colvals))")
    end
    return nothing
end

function b(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    eq_idx, _ = getConstraintTypes(host)
    _, rlb, _ = JuMP.prepProblemBounds(host)
    @assert len == length(eq_idx)
    unsafe_copy!(vec, pointer(rlb[eq_idx]), len)
    println("b ($id) = $(rlb[eq_idx])")
    return nothing
end

function c(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    f, _, _ = JuMP.prepProblemBounds(host)
    @assert len == length(f)
    unsafe_copy!(vec, pointer(f), len)
    println("c ($id) = $(f)")
    return nothing
end

function clow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    _, rlb, _ = JuMP.prepProblemBounds(host)
    _, ineq_idx = getConstraintTypes(host)
    @assert len == length(ineq_idx)
    print("clow ($id) = [")
    for (ind,it) in enumerate(ineq_idx)
        val = (isinf(rlb[it]) ? 0.0 : rlb[it])
        unsafe_store!(vec, convert(Cdouble,val), ind)
        print("$val, ")
    end
    println("]")
    return nothing
end

function cupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    _, ineq_idx = getConstraintTypes(host)
    _, _, rub = JuMP.prepProblemBounds(host)
    @assert len == length(ineq_idx)
    print("cupp ($id) = [")
    for (ind,it) in enumerate(ineq_idx)
        val = (isinf(rub[it]) ? 0.0 : rub[it])
        unsafe_store!(vec, convert(Cdouble,val), ind)
        print("$val, ")
    end
    println("]")
    return nothing
end

function xlow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    @assert len == host.numCols
    print("xlow ($id) = [")
    for it in 1:len
        val = (isinf(host.colLower[it]) ? 0.0 : host.colLower[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
        print("$val, ")
    end
    println("]")
    return nothing
end

function xupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    @assert len == host.numCols
    print("xupp ($id) = [")
    for it in 1:len
        val = (isinf(host.colUpper[it]) ? 0.0 : host.colUpper[it])
        unsafe_store!(vec, convert(Cdouble,val), it)
        print("$val, ")
    end
    println("]")
    return nothing
end

function iclow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    _, ineq_idx = getConstraintTypes(host)
    _, rlb, _ = JuMP.prepProblemBounds(host)
    @assert len == length(ineq_idx)
    print("iclow ($id) = [")
    for (ind,it) in enumerate(ineq_idx)
        val = (isinf(rlb[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), ind)
        print("$val, ")
    end
    println("]")
    return nothing
end

function icupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    _, ineq_idx = getConstraintTypes(host)
    _, _, rub = JuMP.prepProblemBounds(host)
    @assert len == length(ineq_idx)
    print("icupp ($id) = [")
    for (ind,it) in enumerate(ineq_idx)
        val = (isinf(rub[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), ind)
        print("$val, ")
    end
    println("]")
    return nothing
end

function ixlow(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    @assert len == host.numCols
    print("ixlow ($id) = [")
    for it in 1:len
        val = (isinf(host.colLower[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
        print("$val, ")
    end
    println("]")
    return nothing
end

function ixupp(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
    usr = unsafe_pointer_to_objref(user_data)::UserData
    child_id = get_child_index(usr.master, id)
    host = (id == root ? usr.master : usr.children[child_id])
    @assert len == host.numCols
    print("ixupp ($id) = [")
    for it in 1:len
        val = (isinf(host.colUpper[it]) ? 0.0 : 1.0)
        unsafe_store!(vec, convert(Cdouble,val), it)
        print("$val, ")
    end
    println("]")
    return nothing
end

for (mat_name,old_name) in [(:fQ,:Q), (:fA,:A), (:fB,:B), (:fC,:C), (:fD,:D)]
    @eval $mat_name =
            cfunction($old_name, Void, (Ptr{Void},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}))
end

for (nnz_name,old_name) in [(:fnnzQ,:nnzQ), (:fnnzA,:nnzA), (:fnnzB,:nnzB), (:fnnzC,:nnzC), (:fnnzD,:nnzD)]
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
