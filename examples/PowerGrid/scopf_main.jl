include("scopf.jl")

# to obtain a list of off-lines run:  grep "optimal" summaries/case30_SCOPF_summary | awk '{print $2;}' | tr ':\n' ' '


function main()
  if length(ARGS) < 1
    println("Usage: julia scopf_main.jl case_name lines_indexes")
    println("Cases are in 'data' directory: case9 case30 case57 case118 case300 case1354pegase case2383wp case2736sp case2737sop case2746wop case2869pegase case3012wp  case3120sp case3375wp case9241pegase")
    return
  end

  opfdata = opf_loaddata(ARGS[1])

  lines_off=Array(Line, length(ARGS)-1)
  for l in 1:length(lines_off)
    lines_off[l] = opfdata.lines[parse(Int,ARGS[l+1])]
  end
  scopfdata = SCOPFData(opfdata,lines_off)
  @assert length(lines_off) == length(ARGS)-1
  scopfmodel = scopf_model(scopfdata)
  opfmodel,status = scopf_solve(scopfmodel,scopfdata)
  if status==:Optimal
    scopf_outputAll(opfmodel,scopfdata)
  end
end

main()
