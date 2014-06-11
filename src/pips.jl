libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")

function pips_solve(m::Model)
    @assert parent(m) == nothing # make sure this is master problem

    # initialize MPI

    numScens = convert(Cint, length(children(m)))

    function fQ(user_data, id::Cint, krowM::Ptr{Cint}, jcolM::Ptr{Cint}, M::Ptr{Cdouble})
        
    end

    # 
    ccall((libpips,:PIPSSolve), Void, (Ptr{Void},  # MPI_COMM
                                       Cint,       # numScens
                                       Cint,       # nx0
                                       Cint,       # my0
                                       Cint,       # mz0
                                       Cint,       # nx
                                       Cint,       # my
                                       Cint,       # mz
                                       Ptr{Void},  # Q0
                                       Ptr{Void},  # nnzQ0
                                       Ptr{Void},  # c0
                                       Ptr{Void},  # A0
                                       Ptr{Void},  # nnzA0
                                       Ptr{Void},  # B0
                                       Ptr{Void},  # nnzB0
                                       Ptr{Void},  # b0
                                       Ptr{Void},  # C0
                                       Ptr{Void},  # nnzC0
                                       Ptr{Void},  # D0
                                       Ptr{Void},  # nnzD0
                                       Ptr{Void},  # clow0
                                       Ptr{Void},  # iclow0
                                       Ptr{Void},  # cupp0
                                       Ptr{Void},  # icupp0
                                       Ptr{Void},  # xlow0
                                       Ptr{Void},  # ixlow0
                                       Ptr{Void},  # xupp0
                                       Ptr{Void},  # ixupp0
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
                                       fQ0,        
                                       fnnzQ,      
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
