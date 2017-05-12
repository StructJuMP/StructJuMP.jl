#================================================================
 This package solves the following problem parallely

 min  c_0^T x + \sum{i = 1}^S c_i^Ty_i
 s.t. b_0 - A_0 x           \in K_0,
      b_i - A_i x - B_i y_i \in K_i, \forall i = 1,...,S
                x           \in C_0,
                x_I         \in Z 
                        y_i \in C_i, \forall i = 1,...,S

 where input to the Benders engine is:
 c_all is an array of c_0, c_1, c_2, ... , c_S objective coefficients
 A_all is an array of A_0, A_1, A_2, ... , A_S constraint matrices of master variables
 B_all is an array of      B_1, B_2, ... , B_S constraint matrices of scenario variables (note that scenario variables do not appear at master constraints)
 b_all is an array of b_0, b_1, b_2, ... , b_S right hand side of constraints
 K_all is an array of K_0, K_1, K_2, ... , K_S constraint cones
 C_all is an array of C_0, C_1, C_2, ... , C_S variable cones
 v     is an array of master variable types (i.e. :Bin, :Int, :Cont)

 call with: julia -p <num_threads> <myscript>
================================================================#
using JuMP

# this function loads and solves a conic problem and returns its dual
function loadAndSolveConicProblem(c, A, b, K, C, solver)
    
    # load conic model
    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, K, C)
  
    #println(model) 
    #println("process id $(myid()) started")
    #@show c, A, b, K, C
 
    # solve conic model
    MathProgBase.optimize!(model)
    status = MathProgBase.status(model)

    # return status and dual
    #println("process id $(myid()) status $(status)")
    return status, MathProgBase.getdual(model)
end

function loadMasterProblem(c, A, b, K, C, v, num_scen, solver)
    num_var = length(c)
    # load master problem
    master_model = Model(solver=solver)
    @variable(master_model, x[1:num_var])
    for i = 1:num_var
        setcategory(x[i], v[i])
    end
    # add new objective variables for scenarios
    # we assume the objective is sum of nonnegative terms
    @variable(master_model, θ[1:num_scen] >= 0.0)
    @objective(master_model, :Min, sum(c[i]*x[i] for i in 1:num_var) + sum(θ))
    # add constraint cones
    for (cone, ind) in K
        for i in 1:length(ind)
            if cone == :Zero
                @constraint(master_model, dot(A[ind[i],:], x) == b[ind[i]])
            elseif cone == :NonNeg
                @constraint(master_model, dot(A[ind[i],:], x) <= b[ind[i]])
            elseif cone == :NonPos
                @constraint(master_model, dot(A[ind[i],:], x) >= b[ind[i]])
            elseif cone == :SOC
                @constraint(master_model, norm(b[ind[2:end-1]] - A[ind[2:end-1],:]*x) .<= b[ind[1]] - A[ind[1],:]*x)
                @constraint(master_model, dot(A[ind[1],:], x) <= b[ind[1]])
            else
                error("unrecognized cone $cone")
            end
        end
    end       
    # add variable cones
    for (cone, ind) in C
        if cone == :Zero
            for i in ind
                setlowerbound(x[i], 0.0)
                setupperbound(x[i], 0.0)
            end
        elseif cone == :Free
            # do nothing
        elseif cone == :NonNeg
            for i in ind
                setlowerbound(x[i], 0.0)
            end
        elseif cone == :NonPos
            for i in ind
                setupperbound(x[i], 0.0)
            end
        elseif cone == :SOC
            @constraint(master_model, norm(x[ind[2:end-1]]) <= x[ind[1]])
            setlowerbound(x[ind[1]], 0.0)
        else
            error("unrecognized cone $cone")
        end
    end

    #setcategory(master_model, [v;[:Cont for i = 1:num_scen]])
    return master_model, x, θ
end

function addCuttingPlanes(master_model, num_scen, A_all, b_all, output, x, θ, separator, TOL)
    cut_added = false
    # add cutting planes, one per scenario
    for i = 1:num_scen
        #@show typeof(b_all[i+1])
        coef = vec(output[i][2]' * A_all[i+1])
        rhs = vecdot(output[i][2], b_all[i+1])
        #@show size(coef*separator'), size(rhs)
        # add an infeasibility cut
        if output[i][1] == :Infeasible
            # output[i][2] is a ray
            # so alpha * output[i][2] is also valid for any alpha >= 0.
            # Hence output[i][2] might have very large coefficients and alter
            # the numerical accuracy of the master's solver.
            # We scale it to avoid this issue
            scaling = abs(rhs)
            if scaling == 0
              scaling = maximum(coef)
            end
            @constraint(master_model, dot(coef/scaling, x) <= sign(rhs))
            cut_added = true
        # add an optimality cut
        else
            if getvalue(θ[i]) < (dot(coef, separator) - rhs)[1] - TOL
                @constraint(master_model, dot(coef, x) - θ[i] <= rhs)
                cut_added = true
            end
        end
    end
    return cut_added

end

function Benders_pmap(c_all, A_all, B_all, b_all, K_all, C_all, v, master_solver, sub_solver, TOL=1e-5)
    #println("Benders pmap started")
    num_master_var = length(c_all[1])
    num_scen = length(c_all) - 1

    num_bins = zeros(Int, num_scen)
    for i = 1:num_scen
        num_bins[i] = length(b_all[i+1])
    end

    #println("Load master problem")
    (master_model, x, θ) = loadMasterProblem(c_all[1], A_all[1], b_all[1], K_all[1], C_all[1], v, num_scen, master_solver)

    cut_added = true
    separator = zeros(num_master_var)
    objval = Inf
    status = :Infeasible
    while cut_added
        #println("Iteration started")
        status = solve(master_model)
        if status == :Infeasible
            break
        end
        objval = getobjectivevalue(master_model)
        separator = zeros(num_master_var)
        for i = 1:num_master_var
            separator[i] = getvalue(x[i])
        end 
        new_rhs = [zeros(num_bins[i]) for i in 1:num_scen]
        for i = 1:num_scen
            new_rhs[i] = b_all[i+1] - A_all[i+1] * separator
        end

        output = pmap(loadAndSolveConicProblem, 
            [c_all[i+1] for i = 1:num_scen], 
            [B_all[i] for i = 1:num_scen], 
            [new_rhs[i] for i = 1:num_scen], 
            [K_all[i+1] for i = 1:num_scen], 
            [C_all[i+1] for i = 1:num_scen], 
            [sub_solver for i = 1:num_scen])

        #@show output
        cut_added = addCuttingPlanes(master_model, num_scen, A_all, b_all, output, x, θ, separator, TOL)
    end
    return status, objval, separator 
end
