include(string(ENV["HOME"],"/workspace/PIPS/PIPS-NLP/JuliaInterface/ParPipsNlp.jl"))

module ParPipsInterface


using StochJuMP, JuMP
using ParPipsNlp
using MPI

import MathProgBase


export solve


global m
#########################
# Helper 
#########################
function get_model(id)
    return id==0?m:getchildren(m)[id]
end

function get_numvar(id)
    mm = get_model(id)
    nvar_1st = MathProgBase.numvar(m)
    return id==0?nvar_1st:(MathProgBase.numvar(mm)-nvar_1st)
end

function get_numcon(id)
    mm = get_model(id)
    return MathProgBase.numconstr(mm)
end

function array_copy(src,dest)
    @assert length(src)==length(dest)
    array_copy(src,1,dest,1,length(src))
end

function array_copy(src,os, dest, od, n)
    for i=0:n-1
        dest[i+od] = src[i+os]
    end
end

function get_nlp_evaluator(id)
    e = JuMPNLPEvaluator(get_model(id))
    MathProgBase.initialize(e,[:Grad,:Jac,:Hess])
    return e
end

function build_x(id,x0,x1)
    # @show "build_x", id, length(x0), length(x1)
    # @show x0, x1
    if id==0
        return x0
    else
        #build x using index tracking info
        mm = get_model(id)
        othermap = getStochastic(mm).othermap
        new_x = Vector{Float64}(MathProgBase.numvar(mm))
        unsafe_copy!(new_x,1,x1,1,length(x1)) 
        for e in othermap
            pidx = e[1].col
            cidx = e[2].col
            # @show pidx, cidx
            assert(cidx > length(x1))
            new_x[cidx] = x0[pidx]
        end
        return new_x
    end
end

function get_constraints_idx(id)
    lb,ub=getConstraintBounds(get_model(id))
    eq_idx = Dict{Int,Int}()  #jump-> actual used
    ieq_idx = Dict{Int,Int}()
    for i =1:length(lb)
        if lb[i] == ub[i]
            eq_idx[i] = length(eq_idx) + 1
        else
            ieq_idx[i]= length(ieq_idx) + 1
        end
    end
    return (eq_idx,ieq_idx)
end

function get_jac_col_var_idx(rowid, colid)  #this method customerized for no linking constraint presented
    # @show "get_jac_col_var_idx",rowid,colid
    idx_map = Dict{Int,Int}() #dummy (jump) -> actual used
    if rowid == colid
        nvar = get_numvar(rowid)
        for i = 1:nvar
            idx_map[i] = i
        end
    else
        @assert rowid!=0 && rowid != colid
        mm = get_model(rowid)
        othermap = getStochastic(mm).othermap
        for p in othermap
            pidx = p[1].col
            cidx = p[2].col
            idx_map[cidx] = pidx
        end
    end
    # @show idx_map
    return idx_map
end

function convert_to_c_idx(indicies)
    for i in 1:length(indicies)
        indicies[i] = indicies[i] - 1
    end
end

function get_h_var_idx(rowid, colid)
    # @show "get_h_var_idx",rowid,colid
    col_idx_map = Dict{Int,Int}() #dummy (jump) -> actual used
    row_idx_map = Dict{Int,Int}()
    if rowid == colid
        nvar = get_numvar(rowid)
        for i = 1:nvar
            col_idx_map[i] = i
            row_idx_map[i] = i
        end
    elseif rowid == 0  #border
        mm = get_model(colid)
        othermap = getStochastic(mm).othermap
        for p in othermap
            pidx = p[1].col
            cidx = p[2].col
            col_idx_map[cidx] = pidx
        end
        for i = 1:get_numvar(colid)
            row_idx_map[i] = i
        end
    elseif colid == 0 #root contrib.
        mm = get_model(rowid)
        othermap = getStochastic(mm).othermap
        for p in othermap
            pidx = p[1].col
            cidx = p[2].col
            col_idx_map[cidx] = pidx
            row_idx_map[cidx] = pidx
        end
    else
        @assert false rowid colid
    end
    # @show col_idx_map,row_idx_map
    return col_idx_map,row_idx_map
end


function SparseMatrix.sparse(I,J,V, M, N;keepzeros=false)
    if(!keepzeros)
        return sparse(I,J,V,M,N)
    else
        full = sparse(I,J,ones(Float64,length(I)),M,N)
        actual = sparse(I,J,V,M,N)
        fill!(full.nzval,0.0)
        for i in eachindex(actual)
            full[i[1],i[2]] = actual[i]
        end

        return full
    end
end
###############
# Generatic linking code
###############
function str_init_x0(id,x0)
    fill!(x0,1.0)
end

function str_prob_info(id,mode,clb,cub,rlb,rub)
    # @show id
    if mode == :Structure
        nn = get_numvar(id)
        mm = get_numcon(id)
        # @show nn,mm
        return (nn,mm)
    elseif mode == :Values
        # @show length(clb),length(cub)
        mm = get_model(id)
        nvar = get_numvar(id)
        @assert length(clb) == nvar 
        @assert length(cub) == nvar
        array_copy(mm.colUpper, 1, cub, 1, nvar)
        array_copy(mm.colLower, 1, clb, 1, nvar)
        lb,ub = getConstraintBounds(get_model(id))
        (eq_idx, ieq_idx) = get_constraints_idx(id)
        @assert length(lb) == length(ub)
        @assert length(eq_idx) + length(ieq_idx) == length(lb)
        for i in eq_idx
            rlb[i[2]] = lb[i[1]]
            rub[i[2]] = ub[i[1]]
        end
        for i in ieq_idx
            rlb[i[2]] = lb[i[1]]
            rub[i[2]] = ub[i[1]]
        end
    else
        @assert false mode
    end
end

function str_eval_f(id,x0,x1)
    e =  get_nlp_evaluator(id)
    return MathProgBase.eval_f(e,build_x(id,x0,x1))
end

function str_eval_g(id,x0,x1, new_eq_g, new_inq_g)
    e = get_nlp_evaluator(id)
    g = Vector{Float64}(get_numcon(id))
    MathProgBase.eval_g(e,g,build_x(id,x0,x1))
    (eq_idx, ieq_idx) = get_constraints_idx(id)
    @assert length(new_eq_g) == length(eq_idx)
    @assert length(ieq_idx) == length(new_inq_g)
    for i in eq_idx
        new_eq_g[i[2]] = g[i[1]]
    end
    for i in ieq_idx
        new_inq_g[i[2]] = g[i[1]]
    end
end

function str_eval_grad_f(rowid, colid, x0, x1, new_grad_f)
    @assert rowid >= colid
    e = get_nlp_evaluator(rowid)
    x = build_x(rowid,x0,x1)
    g = Vector{Float64}(length(x))
    MathProgBase.eval_grad_f(e,g,x)
    @assert length(g) == MathProgBase.numvar(get_model(rowid))
    @assert length(new_grad_f) == get_numvar(colid)

    var_idx_map = get_jac_col_var_idx(rowid,colid)
    for i in var_idx_map
        new_grad_f[i[2]] = g[i[1]]
    end
end

function str_eval_jac_g(rowid, colid, x0 , x1, mode, e_rowidx, e_colptr, e_values, i_rowidx, i_colptr, i_values)
    # @show "str_eval_jac_g", rowid, colid, mode
    @assert rowid<num_scenarios(m)+1 && colid < num_scenarios(m) + 1
    @assert rowid >= colid
    if(mode == :Structure)
        e = get_nlp_evaluator(rowid)
        (jac_I,jac_J) = MathProgBase.jac_structure(e)
        (eq_idx, ieq_idx) = get_constraints_idx(rowid)
        var_idx = get_jac_col_var_idx(rowid,colid)

        eq_jac_I = Vector{Int}() 
        ieq_jac_I = Vector{Int}()
        eq_jac_J = Vector{Int}()
        ieq_jac_J = Vector{Int}()
        for i = 1:length(jac_I)
            ii = jac_I[i]
            jj = jac_J[i]
            if haskey(var_idx,jj)
                if haskey(eq_idx,ii)
                    push!(eq_jac_I, eq_idx[ii])
                    push!(eq_jac_J, var_idx[jj])
                else
                    @assert haskey(ieq_idx,ii)
                    push!(ieq_jac_I, ieq_idx[ii])
                    push!(ieq_jac_J, var_idx[jj])
                end
            end
        end

        eq_jac = sparse(eq_jac_I,eq_jac_J,ones(Float64,length(eq_jac_I)),length(eq_idx),get_numvar(colid))
        ieq_jac = sparse(ieq_jac_I,ieq_jac_J,ones(Float64,length(ieq_jac_I)),length(ieq_idx),get_numvar(colid))
        # @show (length(eq_jac.nzval), length(ieq_jac.nzval))

        return (length(eq_jac.nzval), length(ieq_jac.nzval))
    elseif(mode == :Values)
        e = get_nlp_evaluator(rowid)
        (jac_I,jac_J) = MathProgBase.jac_structure(e)
        jac_g = Vector{Float64}(length(jac_I))
        MathProgBase.eval_jac_g(e,jac_g,build_x(rowid,x0,x1))
        (eq_idx, ieq_idx) = get_constraints_idx(rowid)
        var_idx = get_jac_col_var_idx(rowid,colid)

        eq_jac_I = Vector{Int}() 
        ieq_jac_I = Vector{Int}()
        eq_jac_J = Vector{Int}()
        ieq_jac_J = Vector{Int}()
        eq_jac_g = Vector{Float64}()
        ieq_jac_g = Vector{Float64}()
        for i = 1:length(jac_I)
            ii = jac_I[i]
            jj = jac_J[i]
            vv = jac_g[i]
            # @show ii, jj, vv
            if haskey(var_idx,jj)
                if haskey(eq_idx,ii)
                    push!(eq_jac_I, eq_idx[ii])
                    push!(eq_jac_J, var_idx[jj])
                    push!(eq_jac_g, vv)
                else
                    @assert haskey(ieq_idx,ii)
                    push!(ieq_jac_I, ieq_idx[ii])
                    push!(ieq_jac_J, var_idx[jj])
                    push!(ieq_jac_g, vv)
                end
            end
        end

        if(length(eq_jac_g) != 0)
            eq_jac = sparse(eq_jac_I,eq_jac_J,eq_jac_g, length(eq_idx),get_numvar(colid), keepzeros=true)
            # @show eq_jac
            # @show eq_jac.rowval
            # @show eq_jac.colptr
            # @show eq_jac.nzval
            array_copy(eq_jac.rowval,1,e_rowidx,1,length(eq_jac.rowval))
            array_copy(eq_jac.colptr,1,e_colptr,1,length(eq_jac.colptr))
            array_copy(eq_jac.nzval, 1,e_values,1,length(eq_jac.nzval))
            convert_to_c_idx(e_rowidx)
            convert_to_c_idx(e_colptr)
        end

        if(length(ieq_jac_g) != 0)
            ieq_jac = sparse(ieq_jac_I,ieq_jac_J,ieq_jac_g, length(ieq_idx),get_numvar(colid), keepzeros=true)
            # @show ieq_jac
            # @show ieq_jac.rowval
            # @show ieq_jac.colptr
            # @show ieq_jac.nzval
            array_copy(ieq_jac.rowval,1,i_rowidx,1,length(ieq_jac.rowval))
            array_copy(ieq_jac.colptr,1,i_colptr,1,length(ieq_jac.colptr))
            array_copy(ieq_jac.nzval, 1,i_values,1,length(ieq_jac.nzval))
            convert_to_c_idx(i_rowidx)
            convert_to_c_idx(i_colptr)
        end
    else
        @assert false mode
    end
    # @show "end str_eval_jac_g"
end


function str_eval_h(rowid, colid, x0, x1, obj_factor, lambda, mode, rowidx, colptr, values)
    # @show "str_eval_h", rowid, colid, mode, length(x0), length(x1),obj_factor, length(lambda), length(rowidx), length(colptr), length(values)
    @assert rowid<num_scenarios(m)+1 && colid < num_scenarios(m) + 1
    low = min(rowid, colid)
    high = max(rowid, colid)
    if(mode == :Structure)
        e = get_nlp_evaluator(high)
        (h_J, h_I) = MathProgBase.hesslag_structure(e) #reverse I, J as we want upper trangular matrix
        col_var_idx, row_var_idx = get_h_var_idx(rowid, colid)
        # @show h_I
        # @show h_J

        new_h_I = Vector{Int}()
        new_h_J = Vector{Int}()
        for i = 1:length(h_I)
            ii = h_I[i]
            jj = h_J[i]
            if haskey(col_var_idx,jj) && haskey(row_var_idx,ii)
                push!(new_h_I,row_var_idx[ii])
                push!(new_h_J,col_var_idx[jj])
            end
        end
        # @show new_h_I
        # @show new_h_J

        laghess = sparse(new_h_I,new_h_J, ones(Float64,length(new_h_I)))
        # @show length(laghess.nzval)
        return length(laghess.nzval)
    elseif(mode == :Values)
        e = get_nlp_evaluator(high)
        (h_J, h_I) = MathProgBase.hesslag_structure(e)
        h = Vector{Float64}(length(h_I))
        MathProgBase.eval_hesslag(e,h,build_x(high,x0,x1),obj_factor,lambda)
        col_var_idx,row_var_idx = get_h_var_idx(rowid, colid)
        
        new_h_I = Vector{Int}()
        new_h_J = Vector{Int}()
        new_h = Vector{Float64}()
        for i = 1:length(h_I)
            ii = h_I[i]
            jj = h_J[i]
            vv = h[i]
            if haskey(col_var_idx,jj) && haskey(row_var_idx,ii)
                push!(new_h_I,row_var_idx[ii])
                push!(new_h_J,col_var_idx[jj])
                push!(new_h, vv)
            end
        end
        # @show new_h_I
        # @show new_h_J
        # @show new_h
        
        if(rowid !=0 && colid == 0)  #root diag contrib.
            # @show "root contrib", rowid, colid
            (h0_I,h0_J) = MathProgBase.hesslag_structure(get_nlp_evaluator(0))
            str_laghess = sparse([new_h_I;h0_I], [new_h_J;h0_J], [new_h;zeros(Float64,length(h0_I))], length(row_var_idx), length(col_var_idx), keepzeros=true)
            # @show str_laghess.m, str_laghess.n, length(str_laghess.nzval)
            # @show str_laghess
            
            array_copy(str_laghess.rowval,1,rowidx,1,length(str_laghess.rowval))
            array_copy(str_laghess.colptr,1,colptr,1,length(str_laghess.colptr))
            array_copy(str_laghess.nzval, 1,values,1,length(str_laghess.nzval))

            # @show str_laghess.nzval
        else
            str_laghess = sparse(new_h_I, new_h_J, new_h, length(row_var_idx), length(col_var_idx), keepzeros=true)
            # @show str_laghess.m, str_laghess.n, length(str_laghess.nzval)
            # @show str_laghess
            
            array_copy(str_laghess.rowval,1,rowidx,1,length(str_laghess.rowval))
            array_copy(str_laghess.colptr,1,colptr,1,length(str_laghess.colptr))
            array_copy(str_laghess.nzval, 1,values,1,length(str_laghess.nzval))
            # @show str_laghess.nzval
        end

        convert_to_c_idx(rowidx)
        convert_to_c_idx(colptr)
    else
        @assert false mode
    end 
    # @show "end str_eval_h"
end


#######
# Linking with PIPS Julia Structure interface
######
function createPipsProblem(model)
    @show "createPipsProblem"

    comm = MPI.COMM_WORLD
    @show "[$(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] create problem "
    @show comm
    
    prob = createProblemStruct(comm,
        num_scenarios(model),  #number scen
        4,  #number var
        3,  #number cons
        str_init_x0, str_prob_info, str_eval_f, str_eval_g, str_eval_grad_f,
        str_eval_jac_g, str_eval_h) 
    @show "end create problem "
    return prob
end

function solvePipsProblem(prob)
    @show "solvePipsProblem"
    @show prob
    ret = solveProblemStruct(prob)
    @show "end solve problem"
end

function solve(model)
    @show "solve"
    global m = model
    MPI.Init()

    prob = createPipsProblem(model)
    
    solvePipsProblem(prob)


    MPI.Finalize()
end

end

