if length(ARGS) == 0 
    Base.warn("No solver is given!")
    @printf "Please provide option for solver to use: Ipopt, ParPips, or Pips \n"
    exit(0);
end
if ARGS[1] == "Ipopt"
    using SerialIpoptInterface
    structJuMPSolve = SerialIpoptInterface.structJuMPSolve
elseif ARGS[1] == "ParPips"
    using ParPipsNlpInterface
    structJuMPSolve = ParPipsNlpInterface.structJuMPSolve
elseif ARGS[1] == "Pips"
    using SerialPipsNlpInterface
    structJuMPSolve = SerialPipsNlpInterface.structJuMPSolve
else
    Base.error("Known solvers are: Ipopt, ParPips, or Pips")
end
