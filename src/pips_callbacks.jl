import Base.Meta.quot

for (mat_name,sym) in [(:fQ,:Q), (:fA,:A), (:fB,:B), (:fC,:C), (:fD,:D)]
    @eval begin
        function $(sym)(user_data::Ptr{Void}, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
            usr = unsafe_pointer_to_objref(user_data)::UserData
            println("id = $id")
            host = (id == root ? usr.master : usr.children[id])
            unsafe_copy!(krowM, pointer(getfield(host,$(quot(sym))).colptr), length(getfield(host,$(quot(sym))).colptr))
            unsafe_copy!(jcolM, pointer(getfield(host,$(quot(sym))).rowval), length(getfield(host,$(quot(sym))).rowval))
            unsafe_copy!(M, pointer(getfield(host,$(quot(sym))).nzval), length(getfield(host,$(quot(sym))).nzval))
            #println($sym, "($id):")
            #println("rowptr    = ", getfield(host,$(quot(sym))).colptr)
            #println("colvals   = ", getfield(host,$(quot(sym))).rowval)
            #println("rownzvals = ", getfield(host,$(quot(sym))).nzval)
            #println()
            return nothing
        end
        $mat_name =
            cfunction($sym, Void, (Ptr{Void},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble}))
    end

end

for (nnz_name,nnzsym,sym) in [(:fnnzQ,:nnzQ,:Q),
                              (:fnnzA,:nnzA,:A),
                              (:fnnzB,:nnzB,:B),
                              (:fnnzC,:nnzC,:C),
                              (:fnnzD,:nnzD,:D)]
    @eval begin
        function $(nnzsym)(user_data::Ptr{Void}, id::Cint, nnz::Ptr{Cint})
            usr = unsafe_pointer_to_objref(user_data)::UserData
            host = (id == root ? usr.master : usr.children[id])
            unsafe_store!(nnz, cint(length(getfield(host,$(quot(sym))).rowval)), 1)
            return nothing
        end
        $nnz_name =
            cfunction($nnzsym, Void, (Ptr{Void},Cint,Ptr{Cint}))
    end
end

for (vec_name,sym) in [(:fb,:b), (:fc,:c),
                            (:fclow,:clow), (:fcupp,:cupp),
                            (:fxlow,:xlow), (:fxupp,:xupp),
                            (:ficlow,:iclow), (:ficupp,:icupp),
                            (:fixlow,:ixlow), (:fixupp,:ixupp)]
    @eval begin
        function $(sym)(user_data::Ptr{Void}, id::Cint, vec::Ptr{Cdouble}, len::Cint)
            usr = unsafe_pointer_to_objref(user_data)::UserData
            host = (id == root ? usr.master : usr.children[id])
            @assert len == length(getfield(host, $(quot(sym))))
            unsafe_copy!(vec, pointer(getfield(host, $(quot(sym)))), len)
            #println($(sym), "($id) = ", getfield(host, $(quot(sym))))
            return nothing
        end
        $vec_name =
            cfunction($sym, Void, (Ptr{Void},Cint,Ptr{Cdouble},Cint))
    end
end
