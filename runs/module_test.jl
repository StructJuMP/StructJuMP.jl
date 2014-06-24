    tic()
    import MPI
    using JuMP
    using StochJuMP
    t1 = toc()

    MPI.init()

    comm = MPI.COMM_WORLD
    root = 0

    MPI.barrier(comm)
    myrank = MPI.rank(comm)

    tt1 = MPI.Reduce(t1, MPI.MAX, root, comm)

    if myrank == 0
        fp = open("$(ENV["HOME"])/.julia/v0.3/StochJuMP/runs/module_load_results.txt", "w")
        println(fp, "tt1 = $tt1 sec")
    end
MPI.finalize()
