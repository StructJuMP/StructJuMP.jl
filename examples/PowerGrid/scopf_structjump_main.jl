include("scopf_structjump.jl")

function main()
  if length(ARGS) < 4
    println("Usage: julia scopf_structjump_main.jl solver enable_profiling case_name lines_indexes")
    println("solver is one of PipsNlp, PipsNlpSerial or Ipopt")
    println("enable_profiling is true or false. If enable_profiling is true and solver is PipsNlp, PIPS should be compiled with NLP_TIMING option.")
    println("Cases are in 'data' directory: case9 case30 case57 case118 case300 case1354pegase case2383wp case2736sp case2737sop case2746wop case2869pegase case3012wp  case3120sp case3375wp case9241pegase")
    return
  end

  solver = ARGS[1]; shift!(ARGS)
  prof = eval(parse(ARGS[1])); shift!(ARGS)
  
  @timing prof tic() #t_user_model_loading
    
  @timing prof tic() #t_rawdata_loading
    
  raw = loadrawdata(ARGS[1])
  @timing prof t_rawdata_loading = toq()
  
  scopfmodel, scopfdata = scopf_model(raw)
  
  # Initial point - needed especially for pegase cases
  scopf_init_x(scopfmodel,scopfdata)
  
  @timing prof t_user_model_loading = toq()
  
  model,status = scopf_solve(scopfmodel,scopfdata, solver, prof)
  if status==:Optimal || status == :Solve_Succeeded
    scopf_structjump_outputAll(scopfmodel,scopfdata)
  end

  @timing prof begin 
    @message @sprintf(" t_user_model_loading %s ", t_user_model_loading)
    @message @sprintf("t_rawdata_loading %s", t_rawdata_loading)   
  end 
end
#run the case 9 example
main()