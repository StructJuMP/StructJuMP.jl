

module SerialIpoptInterface

using StructJuMP, JuMP
using Ipopt

import MathProgBase

include("./structure_helper.jl")

type NonStructJuMPModel
    model::JuMP.Model 
    jac_I::Vector{Int}
    jac_J::Vector{Int}
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    nz_jac::Vector{Int}
    nz_hess::Vector{Int}

    get_x::Function
    numvars::Function
    numcons::Function
    nele_jac::Function
    nele_hess::Function
    get_bounds::Function
    eval_f::Function
    eval_g::Function
    eval_grad_f::Function
    eval_jac_g::Function
    eval_h::Function

    function NonStructJuMPModel(model)
        instance = new(model, 
            Vector{Int}(), Vector{Int}(), Vector{Int}(), Vector{Int}(),
            Vector{Int}(), Vector{Int}()
            )
        
        instance.get_x = function()
            m = instance.model
            v = Float64[];
            @show num_scenarios(m)
            for i = 0:num_scenarios(m)
                mm = get_model(m,i)
                for j = 1:get_numvars(m,i)
                    v_j = getValue(Variable(mm,j))
                    isnan(v_j)? push!(v,1.0):push!(v,v_j)
                end
            end
            @assert length(v) == g_numvars(m)
            @show v
            return v
        end

        instance.numvars = function()
            return g_numvars(instance.model)
        end

        instance.numcons = function()
            return g_numcons(instance.model)
        end

        instance.nele_jac = function()
            @assert length(instance.jac_I) == length(instance.jac_J)
            return length(instance.jac_I)
        end

        instance.nele_hess = function()
            @assert length(instance.hess_I) == length(instance.hess_J)
            return length(instance.hess_I)
        end

        instance.get_bounds = function()
            m = instance.model
            nvar = g_numvars(m)
            ncon = g_numcons(m)
            x_L = Vector{Float64}(nvar)
            x_U = Vector{Float64}(nvar)
            g_L = Vector{Float64}(ncon)
            g_U = Vector{Float64}(ncon)

            row_start = 1
            col_start = 1
            for i = 0:num_scenarios(m)
                mm = get_model(m,i)
                nx = get_numvars(m,i)
                array_copy(mm.colUpper, 1, x_U, col_start, nx)
                array_copy(mm.colLower, 1, x_L, col_start, nx)

                lb,ub = getConstraintBounds(mm)
                ncons = get_numcons(m,i)
                array_copy(lb,1,g_L,row_start,ncons)
                array_copy(ub,1,g_U,row_start,ncons)

                row_start += get_numcons(m,i)
                col_start += get_numvars(m,i)
            end
            return x_L, x_U, g_L, g_U
        end
        
        instance.eval_f = function(x)
            m = instance.model
            objv = 0.0
            start_idx = 1
            for i=0:num_scenarios(m)
                x_new = strip_x(m,i,x,start_idx)    
                objv += MathProgBase.eval_f(get_nlp_evaluator(m,i),x_new)
                start_idx += get_numvars(m,i)
            end
            return objv;
        end 

        instance.eval_g = function(x,g)
            m = instance.model
            @assert length(g) == g_numcons(m)
            start_idx = 1
            g_start_idx = 1
            for i=0:num_scenarios(m)
                x_new = strip_x(instance.model,i,x,start_idx)
                ncon = get_numcons(m,i)    
                g_new = Vector{Float64}(ncon)
                e = get_nlp_evaluator(m,i)
                MathProgBase.eval_g(e,g_new,strip_x(instance.model,i,x,start_idx))
                array_copy(g_new,1,g,g_start_idx,ncon)
                g_start_idx += ncon
                start_idx += get_numvars(m,i)
            end
        end

        instance.eval_grad_f = function(x,grad_f)
            fill!(grad_f,0.0)
            m = instance.model
            start_idx = 1
            for i=0:num_scenarios(m)
                x_new = strip_x(instance.model,i,x,start_idx)
                e = get_nlp_evaluator(m,i)

                g_f = Vector{Float64}(length(x_new))
                MathProgBase.eval_grad_f(e,g_f,x_new)
                nx = get_numvars(m,i) 

                array_copy(g_f,1,grad_f,start_idx,nx)

                othermap = getStructure(get_model(m,i)).othermap
                for i in othermap
                    pid = i[1].col
                    cid = i[2].col
                    grad_f[pid] += g_f[cid]
                    @assert pid <= get_numvars(m,0)
                end
                start_idx += nx
            end
        end

        instance.eval_jac_g = function(x,mode,rows,cols,values)
            m = instance.model
            if mode==:Structure
                @assert length(rows) == length(cols) 
                for i = 1:length(rows)
                    rows[i] = instance.jac_I[i]
                    cols[i] = instance.jac_J[i]
                end
            else
                start_idx = 1
                value_start = 1
                for i = 0:num_scenarios(m)
                    e = get_nlp_evaluator(m,i)
                    x_new = strip_x(instance.model,i,x,start_idx)
                    i_nz_jac = instance.nz_jac[i+1]
                    jac_g = Vector{Float64}(i_nz_jac)
                    MathProgBase.eval_jac_g(e,jac_g,x_new)
                    array_copy(jac_g,1,values,value_start,i_nz_jac)
                    nx = get_numvars(m,i)
                    start_idx += nx
                    value_start += i_nz_jac
                end
            end
        end

        instance.eval_h = function(x, mode, rows, cols, obj_factor, lambda, values)
            m = instance.model
            if mode == :Structure
                @assert length(rows) == length(cols) 
                for i = 1:length(rows)
                    rows[i] = instance.hess_I[i]
                    cols[i] = instance.hess_J[i]
                end
            else
                start_idx = 1
                value_start = 1
                lambda_start = 1
                for i = 0:num_scenarios(m)
                    e = get_nlp_evaluator(m,i)
                    x_new = strip_x(instance.model,i,x,start_idx)
                    nc = get_numcons(m,i)
                    lambda_new = Vector{Float64}(nc)
                    array_copy(lambda,lambda_start, lambda_new, 1, nc)
                    i_nz_hess = instance.nz_hess[i+1]
                    h = Vector{Float64}(i_nz_hess)
                    MathProgBase.eval_hesslag(e,h,x_new,obj_factor,lambda_new)
                    array_copy(h,1,values,value_start,i_nz_hess)
                    nx = get_numvars(m,i)
                    start_idx += nx
                    lambda_start += nc
                    value_start += i_nz_hess
                end
            end
        end

        # initialization  jac
        col_offset = 0
        row_offset = 0
        m = instance.model
        for i = 0:num_scenarios(m)
            reverse_map = Dict{Int,Int}()
            mm = get_model(m,i)
            for ety in getStructure(mm).othermap
                reverse_map[ety[2].col] = ety[1].col #child->parent
            end
            # @show reverse_map
            e = get_nlp_evaluator(m,i)
            # @show "after e"
            i_jac_I, i_jac_J =  MathProgBase.jac_structure(e)
            # @show "after strct jac"
            for idx = 1:length(i_jac_J)
                jj = i_jac_J[idx]
                if haskey(reverse_map,jj)
                    push!(instance.jac_J, reverse_map[jj])
                else
                    push!(instance.jac_J, jj + col_offset)
                end
                push!(instance.jac_I, i_jac_I[idx] + row_offset)
            end
            push!(instance.nz_jac, length(i_jac_J)) #offset by 1

            col_offset += get_numvars(m,i)
            row_offset += get_numcons(m,i)
        end

        #initialization hess
        offset = 0
        for i = 0:num_scenarios(m)
            reverse_map = Dict{Int,Int}()
            mm = get_model(m,i)
            for ety in getStructure(mm).othermap
                reverse_map[ety[2].col] = ety[1].col #child->parent
            end

            e = get_nlp_evaluator(m,i)
            i_hess_I, i_hess_J =  MathProgBase.hesslag_structure(e)
            for idx = 1:length(i_hess_I)
                ii = i_hess_I[idx]
                jj = i_hess_J[idx]

                if haskey(reverse_map,ii)
                    push!(instance.hess_I, reverse_map[ii])
                else
                    push!(instance.hess_I, ii + offset)
                end

                if haskey(reverse_map,jj)
                    push!(instance.hess_J, reverse_map[jj])
                else
                    push!(instance.hess_J, jj + offset)
                end
            end
            push!(instance.nz_hess, length(i_hess_I)) #offset by 1

            offset += get_numvars(m,i)
        end

        return instance  
    end
end

function structJuMPSolve(model; suppress_warmings=false,kwargs...)
    # @show typeof(model)
    nm = NonStructJuMPModel(model)
    x_L, x_U, g_L, g_U = nm.get_bounds()
    n = g_numvars(model)
    m = g_numcons(model)
    nele_jac = nm.nele_jac()
    nele_hess = nm.nele_hess()

    # @show x_L, x_U
    # @show g_L, g_U
    # @show n,m
    # @show nele_jac,nele_hess

    prob = createProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                         nm.eval_f, nm.eval_g, nm.eval_grad_f, nm.eval_jac_g, nm.eval_h)

    prob.x = nm.get_x()
    status = solveProblem(prob)
    return status
end

end