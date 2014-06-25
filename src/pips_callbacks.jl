for (sym,nnzsym) in [(:Q,:nnzQ), (:A,:nnzA), (:B,:nnzB), (:C,:nnzC), (:D,:nnzD)]
    @eval begin
        function $(sym)(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
            usr = unsafe_pointer_to_objref(user_data)::UserData
            host = (id == root ? usr.master : usr.children[get_child_index(usr.master, id)])
            unsafe_copy!(krowM, host.$(sym).rowptr.-1, length(host.$(sym).rowptr))
            unsafe_copy!(jcolM, host.$(sym).colvals.-1, length(host.$(sym).colvals))
            unsafe_copy!(M, host.$(sym).nzvals, length(host.$(sym).nzvals))
            println("$(sym) ($id):")
            println("rowptr    = $($(sym).rowptr)")
            println("colvals   = $($(sym).colvals)")
            println("rownzvals = $($(sym).nzvals)")
            println()
            return nothing            
        end
        function $(nnzsym)(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
            usr = unsafe_pointer_to_objref(user_data)::UserData
            host = (id == root ? usr.master : usr.children[get_child_index(usr.master, id)])
            unsafe_store!(nnz, cint(length(host.$(sym).colvals)), 1)
            return nothing
        end
    end
end

for (sym) in [:b,:c,:clow,:cupp,:xlow,:xupp,:iclow,:icupp,:ixlow,:ixupp]
    @eval begin
        function $(sym)(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
            usr = unsafe_pointer_to_objref(user_data)::UserData
            host = (id == root ? usr.master : usr.children[get_child_index(usr.master, id)])
            @assert len == length($(sym))
            unsafe_copy!(vec, pointer(host.$(sym)), len)
            println("$(sym) ($id) = $(host.$(sym))")
            return nothing
    end
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
