include("opfdata_structjump.jl")
using StructJuMP, JuMP
using StructJuMPSolverInterface

type SCOPFData
  raw::RawData
  lines_off::Array
  #Float64::gener_ramp #generator ramp limit for contingency (percentage)
end

function scopf_solve(scopfmodel, scopfdata, solver, prof::Bool) 
  
  status = solve(scopfmodel;solver=solver, with_prof=prof)
  
  # getVarValue(scopfmodel)
  return scopfmodel,status
end

function scopf_model(raw::RawData)
  opfdata = opf_loaddata(raw) #load root node
  lines_off=Array(Line, length(ARGS)-1)
  for l in 1:length(lines_off)
    lines_off[l] = opfdata.lines[parse(Int,ARGS[l+1])]
  end
  sd = SCOPFData(raw,lines_off)

  nscen = length(sd.lines_off)

  opfmodel = StructuredModel(num_scenarios=nscen)

  #shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA;
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus  = length(buses); nline = length(lines); ngen  = length(generators);
 
  #branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)

  @variable(opfmodel, generators[i].Pmin <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel, 0<=extra[i=1:ngen]<=0.05*generators[i].Pmax)
  @variable(opfmodel, generators[i].Qmin <= Qg[i=1:ngen] <= generators[i].Qmax)

  @variable(opfmodel, buses[i].Vmin <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel,Va[1:nbus])
  #fix the voltage angle at the reference bus
  setlowerbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
  setupperbound(Va[opfdata.bus_ref], buses[opfdata.bus_ref].Va)
  # setlowerbound(Va[opfdata.bus_ref], -.1)
  # setupperbound(Va[opfdata.bus_ref], .1)

  #objective function
  @NLobjective(opfmodel, Min, (1/(nscen+1))*sum{ generators[i].coeff[generators[i].n-2]*(baseMVA*(Pg[i] + extra[i]))^2 
                         +generators[i].coeff[generators[i].n-1]*(baseMVA*(Pg[i]+extra[i]))
                     +generators[i].coeff[generators[i].n  ], i=1:ngen})


  # @show "Root - Buses: %d  Lines: %d  Generators", nbus, nline, ngen
  
  # power flow balance
  for b in 1:nbus
    #real part
    @NLconstraint(
      opfmodel, 
      ( sum{ YffR[l], l in FromLines[b]} + sum{ YttR[l], l in ToLines[b]} + YshR[b] ) * Vm[b]^2 
      + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }  
      + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   } 
      - ( sum{baseMVA*(Pg[g]+extra[g]), g in BusGeners[b]} - buses[b].Pd ) / baseMVA      # Sbus part
      ==0)

    #imaginary part
    @NLconstraint(
      opfmodel,
      ( sum{-YffI[l], l in FromLines[b]} + sum{-YttI[l], l in ToLines[b]} - YshI[b] ) * Vm[b]^2 
      + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }
      + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   }
      - ( sum{baseMVA*Qg[g], g in BusGeners[b]} - buses[b].Qd ) / baseMVA      #Sbus part
      ==0)
  end

  # generator max
  @NLconstraint(opfmodel, ex[i=1:ngen], generators[i].Pmin <= Pg[i] + extra[i] <= generators[i].Pmax)

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
      @NLconstraint(opfmodel,
        Vm[busIdx[lines[l].from]]^2 *
        ( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2 
          + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])) 
        ) 
        - flowmax <=0)
      
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
  
  # @show "Lines with limits  ", nlinelim

  for c in getLocalChildrenIds(opfmodel) 
  # @declare_second_stage(opfmodel, c, begin
    opfmodel_c = StructuredModel(parent=opfmodel,id=c)
    opfdata_c=opf_loaddata(raw, sd.lines_off[c]) 
    
    #shortcuts for compactness
    lines = opfdata_c.lines
    buses = opfdata_c.buses
    generators = opfdata_c.generators 
    baseMVA = opfdata_c.baseMVA

    busIdx = opfdata_c.BusIdx
    FromLines = opfdata_c.FromLines
    ToLines = opfdata_c.ToLines
    BusGeners = opfdata_c.BusGenerators

    nbus  = length(buses)
    nline = length(lines)
    ngen  = length(generators)
  
    # @printf("Contingency %d -> Buses: %d  Lines: %d  Generators: %d\n", c, nbus, nline, ngen)

    YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(opfdata_c.lines, opfdata_c.buses, opfdata_c.baseMVA)

    @variable(opfmodel_c, 0<=extra[i=1:ngen]<=0.05*generators[i].Pmax)
    @variable(opfmodel_c, buses[i].Vmin <= Vm[i=1:nbus] <= buses[i].Vmax)
    @variable(opfmodel_c, Va[1:nbus])

    #fix the voltage angle at the reference bus
    setlowerbound(Va[opfdata_c.bus_ref], buses[opfdata_c.bus_ref].Va)
    setupperbound(Va[opfdata_c.bus_ref], buses[opfdata_c.bus_ref].Va)

    @NLobjective(opfmodel_c, Min, (1/(nscen+1))*sum{ generators[i].coeff[generators[i].n-2]*(baseMVA*(Pg[i] + extra[i]))^2 
                         +generators[i].coeff[generators[i].n-1]*(baseMVA*(Pg[i]+extra[i]))
                     +generators[i].coeff[generators[i].n  ], i=1:ngen})

    
    # power flow balance
    for b in 1:nbus
      #real part
      @NLconstraint(
        opfmodel_c, 
        ( sum{ YffR[l], l in FromLines[b]} + sum{ YttR[l], l in ToLines[b]} + YshR[b] ) * Vm[b]^2 
        + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }  
        + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   } 
        - ( sum{baseMVA*(Pg[g] + extra[g]), g in BusGeners[b]} - buses[b].Pd ) / baseMVA      # Sbus part
        ==0)

      #imaginary part
      @NLconstraint(
        opfmodel_c,
        ( sum{-YffI[l], l in FromLines[b]} + sum{-YttI[l], l in ToLines[b]} - YshI[b] ) * Vm[b]^2 
        + sum{ Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )), l in FromLines[b] }
        + sum{ Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])), l in ToLines[b]   }
        - ( sum{baseMVA*Qg[g], g in BusGeners[b]} - buses[b].Qd ) / baseMVA      #Sbus part
        ==0)
    end

    #generator max
    @NLconstraint(opfmodel_c, ex[i=1:ngen], generators[i].Pmin <= Pg[i] + extra[i] <= generators[i].Pmax)

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
          opfmodel_c,
          Vm[busIdx[lines[l].from]]^2 *
          ( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2 
            + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])) 
          ) 
          - flowmax <=0)
  
        #branch apparent power limits (to bus)
        Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
        Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
        @NLconstraint(
          opfmodel_c, 
          Vm[busIdx[lines[l].to]]^2 *
          ( Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
            + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
          )
          - flowmax <=0)
      end
    end

    # @printf("Contingency %d -> Buses: %d  Lines: %d  Generators: %d\n", c, nbus, nline, ngen)
    # println("     lines with limits:  ", nlinelim)
  end
  # ) #end second stage
  # @show "acopf_model_sc  - done"
  return opfmodel, sd
end


function scopf_structjump_outputAll(opfmodel, scopf_data)
  #shortcuts for compactness
  sd=scopf_data; opf_data = opf_loaddata(sd.raw);
  lines = opf_data.lines; buses = opf_data.buses; generators = opf_data.generators; baseMVA = opf_data.baseMVA
  busIdx = opf_data.BusIdx; FromLines = opf_data.FromLines; ToLines = opf_data.ToLines; BusGeners = opf_data.BusGenerators;

  nbus  = length(buses); nline = length(lines); ngen  = length(generators)

  # OUTPUTING
  println("Objective value: ", getObjectiveVal(opfmodel), "USD/hr")
  VM=getvalue(getvariable(opfmodel,:Vm)); VA=getvalue(getvariable(opfmodel,:Va))
  PG=getvalue(getvariable(opfmodel,:Pg)); QG=getvalue(getvariable(opfmodel,:Qg))

  EX=getvalue(getvariable(opfmodel,:extra));

  # printing the first stage variables
  println("============================= BUSES ==================================")
  println("  Generator  |  extra ")   # |    P (MW)     Q (MVAr)")  #|         (load)   ")  
  println("----------------------------------------------------------------------")
  for i in 1:ngen
      @printf("  %10d | %6.2f \n",generators[i].bus, EX[i])
  end
  println("\n")


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
  return 

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

  return
end


function scopf_init_x(scopfmodel,scopfdata)
  raw = scopfdata.raw
  lines_off = scopfdata.lines_off
  for i in getLocalBlocksIds(scopfmodel)
    if(i==0)
      opfdata = opf_loaddata(raw)
      Pg0,Qg0,Vm0,Va0 = scopf_compute_x0(opfdata)
      extra0 = 0.025*Pg0
      setvalue(getvariable(scopfmodel, :Pg), Pg0)
      setvalue(getvariable(scopfmodel, :extra), extra0)  
      setvalue(getvariable(scopfmodel, :Qg), Qg0)
      setvalue(getvariable(scopfmodel, :Vm), Vm0)
      setvalue(getvariable(scopfmodel, :Va), Va0)
    else
      mm = getchildren(scopfmodel)[i]
      opfdata_c=opf_loaddata(raw,lines_off[i]) 
      Pg0,Qg0,Vm0,Va0 = scopf_compute_x0(opfdata_c)
      extra0 = 0.025*Pg0
      setvalue(getvariable(mm, :extra), extra0)
      setvalue(getvariable(mm, :Vm), Vm0)
      setvalue(getvariable(mm, :Va), Va0)
    end
  end
end


# Compute initial point for IPOPT based on the values provided in the case data
function scopf_compute_x0(opfdata)
  Pg=zeros(length(opfdata.generators)); 
  Qg=zeros(length(opfdata.generators)); 
  i=1
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
