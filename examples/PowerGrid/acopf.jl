include("opfdata.jl")
using JuMP
using Ipopt

function acopf_solve(opfmodel, opf_data)
   
  # 
  # Initial point - needed especially for pegase cases
  #
  Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opf_data)
  setvalue(getvariable(opfmodel, :Pg), Pg0)  
  setvalue(getvariable(opfmodel, :Qg), Qg0)
  setvalue(getvariable(opfmodel, :Vm), Vm0)
  setvalue(getvariable(opfmodel, :Va), Va0)

  status = solve(opfmodel)

  if status != :Optimal
    println("Could not solve the model to optimality.")
  end
  return opfmodel,status
end

function acopf_model(opf_data)
  #shortcuts for compactness
  lines = opf_data.lines; buses = opf_data.buses; generators = opf_data.generators; baseMVA = opf_data.baseMVA
  busIdx = opf_data.BusIdx; FromLines = opf_data.FromLines; ToLines = opf_data.ToLines; BusGeners = opf_data.BusGenerators;

  nbus  = length(buses); nline = length(lines); ngen  = length(generators)

  #branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)

  #
  # JuMP model now
  #
  opfmodel = Model(solver=IpoptSolver(print_level=0))

  @variable(opfmodel, generators[i].Pmin <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel, generators[i].Qmin <= Qg[i=1:ngen] <= generators[i].Qmax)

  @variable(opfmodel, buses[i].Vmin <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel, Va[1:nbus])
  #fix the voltage angle at the reference bus
  setlowerbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
  setupperbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)

  # minimize active power
#  @NLobjective(opfmodel, 
#		  Min, 
#		  sum{ generators[i].coeff[generators[i].n] + 
#		       sum{generators[i].coeff[generators[i].n-k]*(baseMVA*Pg[i])^k, k=1:generators[i].n-1}, 
#		       i=1:ngen}
#		 )
 
  @NLobjective(opfmodel, Min, sum{ generators[i].coeff[generators[i].n-2]*(baseMVA*Pg[i])^2 
			             +generators[i].coeff[generators[i].n-1]*(baseMVA*Pg[i])
				     +generators[i].coeff[generators[i].n  ], i=1:ngen})

  #
  # power flow balance
  #
  for b in 1:nbus
    #real part
    @NLconstraint(
      opfmodel, 
      ( sum{ YffR[l], l in FromLines[b]} + sum{ YttR[l], l in ToLines[b]} + YshR[b] ) * Vm[b]^2 
      + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }  
      + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   } 
      - ( sum{baseMVA*Pg[g], g in BusGeners[b]} - buses[b].Pd ) / baseMVA      # Sbus part
      ==0)
#  end
#  for b in 1:nbus 
    #imaginary part
    @NLconstraint(
      opfmodel,
      ( sum{-YffI[l], l in FromLines[b]} + sum{-YttI[l], l in ToLines[b]} - YshI[b] ) * Vm[b]^2 
      + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }
      + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   }
      - ( sum{baseMVA*Qg[g], g in BusGeners[b]} - buses[b].Qd ) / baseMVA      #Sbus part
      ==0)
  end
  #
  # branch/lines flow limits
  #
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      #branch apparent power limits (from bus)
      Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
      Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
      @NLconstraint(
        opfmodel,
	Vm[busIdx[lines[l].from]]^2 *
	( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2 
	  + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])) 
	) 
        - flowmax <=0)
#    end
#  end
#  for l in 1:nline
#    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
#      
#      flowmax=(lines[l].rateA/baseMVA)^2

      #branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      @NLconstraint(
        opfmodel,
	Vm[busIdx[lines[l].to]]^2 *
        ( Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
          + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
        )
        - flowmax <=0)
    end
  end
  
  @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
  println("Lines with limits  ", nlinelim)
 
  return opfmodel
end

  #######################################################
  
  #values = zeros(2*nbus+2*ngen) 
  ## values[1:2*nbus+2*ngen] = readdlm("/sandbox/petra/work/installs/matpower5.1/vars2.txt")
  #values[1:2*nbus+2*ngen] = readdlm("/sandbox/petra/work/installs/matpower5.1/vars3_9241.txt")
  #d = JuMP.NLPEvaluator(opfmodel)
  #MathProgBase.initialize(d, [:Jac])

  #g = zeros(2*nbus+2*nlinelim)
  #MathProgBase.eval_g(d, g, values)
  #println("f=", MathProgBase.eval_f(d,values))
 
  #gmat=zeros(2*nbus+2*nlinelim)
  #gmat[1:end] = readdlm("/sandbox/petra/work/installs/matpower5.1/cons3_9241.txt")
  #println("diff: ", norm(gmat-g))

  #println(opfmodel)

  #############################################################

function acopf_outputAll(opfmodel, opf_data)
  #shortcuts for compactness
  lines = opf_data.lines; buses = opf_data.buses; generators = opf_data.generators; baseMVA = opf_data.baseMVA
  busIdx = opf_data.BusIdx; FromLines = opf_data.FromLines; ToLines = opf_data.ToLines; BusGeners = opf_data.BusGenerators;

  nbus  = length(buses); nline = length(lines); ngen  = length(generators)

  # OUTPUTING
  println("Objective value: ", getobjectivevalue(opfmodel), "USD/hr")
  VM=getvalue(getvariable(opfmodel,:Vm)); VA=getvalue(getvariable(opfmodel,:Va))
  PG=getvalue(getvariable(opfmodel,:Pg)); QG=getvalue(getvariable(opfmodel,:Qg))

  println("============================= BUSES ==================================")
  println("  BUS    Vm     Va   |   Pg (MW)    Qg(MVAr) ")   # |    P (MW)     Q (MVAr)")  #|         (load)   ") 
  
  println("                     |     (generation)      ") 
  println("----------------------------------------------------------------------")
  for i in 1:nbus
    @printf("%4d | %6.2f  %6.2f | %s  | \n",
	    buses[i].bus_i, VM[i], VA[i]*180/pi, 
	    length(BusGeners[i])==0?"   --          --  ":@sprintf("%7.2f     %7.2f", baseMVA*PG[BusGeners[i][1]], baseMVA*QG[BusGeners[i][1]]))
  end   
  println("\n")

  within=20 # percentage close to the limits
  
  
  nflowlim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nflowlim += 1
    end
  end

  if nflowlim>0 
    println("Number of lines with flow limits: ", nflowlim)

    optvec=zeros(2*nbus+2*ngen)
    optvec[1:ngen]=PG
    optvec[ngen+1:2*ngen]=QG
    optvec[2*ngen+1:2*ngen+nbus]=VM
    optvec[2*ngen+nbus+1:2*ngen+2*nbus]=VA

    d = JuMP.NLPEvaluator(opfmodel)
    MathProgBase.initialize(d, [:Jac])

    consRhs = zeros(2*nbus+2*nflowlim)
    MathProgBase.eval_g(d, consRhs, optvec)  


    #println(consRhs)

    @printf("================ Lines within %d %s of flow capacity ===================\n", within, "\%")
    println("Line   From Bus    To Bus    At capacity")

    nlim=1
    for l in 1:nline
      if lines[l].rateA!=0 && lines[l].rateA<1.0e10
        flowmax=(lines[l].rateA/baseMVA)^2
        idx = 2*nbus+nlim
        
        if( (consRhs[idx]+flowmax)  >= (1-within/100)^2*flowmax )
          @printf("%3d      %3d      %3d        %5.3f%s\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx]+flowmax)/flowmax), "\%" ) 
          #@printf("%7.4f   %7.4f    %7.4f \n", consRhs[idx], consRhs[idx]+flowmax,  flowmax)
        end
        nlim += 1
      end
    end
  end

  #println(getvalue(Vm))
  #println(getvalue(Va)*180/pi)

  #println(getvalue(Pg))
  #println(getvalue(Qg))

  return
end


# Compute initial point for IPOPT based on the values provided in the case data
function acopf_initialPt_IPOPT(opfdata)
  Pg=zeros(length(opfdata.generators)); Qg=zeros(length(opfdata.generators)); i=1
  for g in opfdata.generators
    # set the power levels in in between the bounds as suggested by matpower 
    # (case data also contains initial values in .Pg and .Qg - not used with IPOPT)
    Pg[i]=0.5*(g.Pmax+g.Pmin)
    Qg[i]=0.5*(g.Qmax+g.Qmin)
    i=i+1
  end
  @assert i-1==length(opfdata.generators)

  Vm=zeros(length(opfdata.buses)); i=1;
  for b in opfdata.buses
    # set the ini val for voltage magnitude in between the bounds 
    # (case data contains initials values in Vm - not used with IPOPT)
    Vm[i]=0.5*(b.Vmax+b.Vmin); 
    i=i+1
  end
  @assert i-1==length(opfdata.buses)

  # set all angles to the angle of the reference bus
  Va = opfdata.buses[opfdata.bus_ref].Va * ones(length(opfdata.buses))

  return Pg,Qg,Vm,Va
end

