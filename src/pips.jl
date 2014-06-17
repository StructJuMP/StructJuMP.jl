libpips = dlopen("/home/huchette/PIPS/PIPS/build/PIPS-IPM/libpipsipm-shared.so")
PIPSSolve = dlsym(libpips,:PIPSSolve)

type UserData
    master   :: JuMP.Model
    children :: Vector{JuMP.Model}
end

const root = 0

cint(a::Int) = convert(Cint,a)
vcint(a::Vector{Int}) = pointer(convert(Vector{Cint},a))::Ptr{Cint}

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

    return rowptr, colval, rownzval
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
    Q = sparse(vars1, vars2, coeff_copy, n, n)
    istriu(Q) && (Q = Q')
    @assert istril(Q)
    return Q.colptr, Q.rowval, Q.nzval
end

include("pips_callbacks.jl")

function pips_solve(master::JuMP.Model)
    @assert getparent(master) == nothing # make sure this is master problem
    @assert master.objSense == :Min

    children = getchildren(master)
    child    = children[1]

    # MPI data
    comm = MPI.COMM_WORLD

    size = MPI.size(comm)
    rank = MPI.rank(comm)

    numScens = num_scenarios(master)
    scenPerRank = iceil(numScens/size)
    println("length(children) = $(length(children))")
    println("num scenarios    = $numScens")
    println("scens per rank   = $scenPerRank")
    # @assert (length(children) == numScens == scenPerRank) # while we're running on one proc

    user_data = UserData(master, children)

    eq_idx_m, ineq_idx_m = getConstraintTypes(master)
    eq_idx_c, ineq_idx_c = getConstraintTypes(child) # should probably check this is const across children

    n_eq_m, n_ineq_m = length(eq_idx_m), length(ineq_idx_m)
    n_eq_c, n_ineq_c = length(eq_idx_c), length(ineq_idx_c)

    mpi_comm = Cint[comm.fval]

    MPI.barrier(comm)

    println("comm (julia) = $(comm.fval)")

    obj_val = [0.0]
    first_primal  = Array(Cdouble, master.numCols)
    second_primal = Array(Cdouble, numScens*child.numCols)
    first_dual    = Array(Cdouble, n_ineq_m)
    second_dual   = Array(Cdouble, numScens*n_eq_m)

    val = ccall(PIPSSolve, Void, (Ptr{Cint},  # MPI_COMM
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

    println("objective = $obj_val")
    println("first stage primal sol  = $first_primal")
    println("second stage primal sol = $second_primal")
    println("first stage dual sol    = $first_dual")
    #println("second stage dual sol   = $second_dual")

    MPI.finalize()
end
