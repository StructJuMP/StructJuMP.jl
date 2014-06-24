function foo()
    2 + 3
    tic()
    import MPI
    t1 = toc()
    tic()
    using JuMP
    t2 = toc()
    tic()
    using JuMPStoch
    t3 = toc()

    comm = MPI.COMM_WORLD
    myrank = MPI.rank(comm)

    if myrank == 0
        println("t1 = $t1 sec")
        println("t2 = $t2 sec")
        println("t3 = $t3 sec")
    end
end
