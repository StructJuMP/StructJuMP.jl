include("acopf.jl")

function main()
  if length(ARGS) != 1
    println("Usage: julia acopf.jl case_name")
    println("Cases are in 'data' directory: case9 case30 case57 case118 case300 case1354pegase case2383wp case2736sp case2737sop case2746wop case2869pegase case3012wp  case3120sp case3375wp case9241pegase")
    return
  end

  opfdata = opf_loaddata(ARGS[1])
  opfmodel = acopf_model(opfdata)
  opfmodel,status = acopf_solve(opfmodel,opfdata)
  if status==:Optimal
    acopf_outputAll(opfmodel,opfdata)
  end
end

main()
