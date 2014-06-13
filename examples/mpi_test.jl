using StochJuMP

MPI.init()
comm = MPI.COMM_WORLD

println("rank = $(MPI.rank(comm))")

MPI.finalize()
