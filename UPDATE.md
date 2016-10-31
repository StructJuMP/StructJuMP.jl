# 08-03-2016

- Removed unnecessary dependencies to ``StructJuMPSolverInterface``, ``PIPS-NLP`` and ``Ipopt``. The current version of ``StructJuMPSolverInterface`` is specific to PIPS and Ipopt. (needs to be generalized)
- Removed requirement of ``MPI.jl``. This makes the ``StructJuMP`` package general enough for any user.
- Removed dependency to ``JuMP`` by exporting all JuMP functions. So, user does not need to say ``using JuMP``.
- Updated README.md. All solver specific statements are removed. They should be moved to each solver interface. Instead, three solvers that can read models from ``StructJuMP`` are briefly introduced.