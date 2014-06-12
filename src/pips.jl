import MPI

libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")

function pips_solve(m::Model)
    @assert parent(m) == nothing # make sure this is master problem

    # initialize MPI
    MPI.init()
    comm = MPI.COMM_WORLD

    numScens = num_scenarios(m)

    f, rowlb, rowub = JuMP.prepProblemBounds(m)

    eq_idx, ineq_idx = getConstraintTypes(m)
    n_eq, n_ineq = length(eq_idx), length(ineq_idx)
    rlb, rub = rowlb[eq_idx], rowub[ineq_idx]

    #####################################################
    # Callback functions for matrices, vectors, and nnz's
    #####################################################
    function fQ(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        unsafe_copy!(krowM, Q_rowptr[id], N)
    end

    for (name, typ) in [(:fA, "eq"),
                        (:fB, "eq"), 
                        (:fC, "ineq"), 
                        (:fD, "ineq")]
        @eval begin
            function $(name)(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
                unsafe_copy!(krowM, symbol(typ*"_rowptr[id]"),    n_eq+1)
                unsafe_copy!(jcolM, symbol(typ*"_colvals[id]"),   symbol("length("*typ*"_colvals[id])"))
                unsafe_copy!(M,     symbol(typ*"_rownzvals[id]"), symbol("length("*typ*"_colvals[id])"))
                return C_NULL
            end
        end
    end

    for (name, lngth) in [(:fnnzQ, :(length(m.obj.qvars1))), 
                          (:fnnzA, :(length(eq_colvals[id]))),
                          (:fnnzB, :(length(eq_colvals[id]))), 
                          (:fnnzC, :(length(eq_colvals[id]))),
                          (:fnnzD, :(length(eq_colvals[id])))]
        @eval begin
            function $(name)(user_data, id::Cint, nnz::Ptr{Cint})
                unsafe_store!(nnz, $(lngth), 1)
                return C_NULL
        end
    end

    for (name, lngth, stor) in [(:fb, :n_eq, :eq_rhs),
                                (:fc, :(length(f)), :f)]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                unsafe_copy!(vec, $(stor), len)
                return C_NULL
            end
        end
    end

    for (name, lngth, stor) in [(:fclow, :n_ineq, :rlb),
                                (:fcupp, :n_ineq, :rub),
                                (:fxlow, :(m.numCols), :(m.colLower)),
                                (:fxupp, :(m.numCols), :(m.colUpper))]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                for it in 1:len
                    val = (isinf($stor)[it]) ? 0.0 : $(stor)[it]
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
            end
        end
    end

    for (name, lngth, stor) in [(:ficlow, :n_ineq, :rlb),
                                (:ficupp, :n_ineq, :rub),
                                (:fixlow, :(m.numCols), :(m.colLower)),
                                (:fixupp, :(m.numCols), :(m.colUpper))]
        @eval begin
            function $(name)(user_data, id::Cint, vec::Ptr{Cdouble}, len::Cint)
                @assert len == $(lngth)
                for it in 1:len
                    val = (isinf($(stor)[it]) ? 0.0 : 1.0)
                    unsafe_store!(vec, convert(Cdouble,val), it)
                end
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
                                      (comm,       
                                       numScens,   
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
