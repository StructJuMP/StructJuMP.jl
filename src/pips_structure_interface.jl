try
    include(ENV["PIPS_NLP_PAR_JULIA_INTERFACE"])
catch err
    if(isa(err, ErrorException))
      warn("Could not include PIPS-NLP Julia interface file. Please setup ENV variable 'PIPS_NLP_PAR_JULIA_INTERFACE' to the location of this file, usually in PIPS repo at PIPS-NLP/JuliaInterface/ParPipsNlp.jl")
    end
    rethrow()
end

module ParPipsInterface


using StructJuMP, JuMP
using ParPipsNlp
using MPI

import MathProgBase

export solve

type StructJuMPModel <: ModelInterface
    internalModel::JuMP.Model
    status::Int
    id_con_idx_map::Dict{Int,Pair{Dict{Int,Int},Dict{Int,Int}}}  #eq, ieq, jump->actual used
    
    get_num_scen::Function
    get_sense::Function
    get_status::Function

    get_num_rows::Function
    get_num_cols::Function
    get_num_eq_cons::Function
    get_num_ineq_cons::Function

    set_status::Function
    
    set_num_rows::Function
    set_num_cols::Function
    set_num_eq_cons::Function
    set_num_ineq_cons::Function

    str_init_x0::Function
    str_prob_info::Function
    str_eval_f::Function
    str_eval_g::Function
    str_eval_grad_f::Function
    str_eval_jac_g::Function
    str_eval_h::Function

    function StructJuMPModel(model::JuMP.Model, status::Int)
        instance = new(model,status,Dict{Int,Pair{Dict{Int,Int},Dict{Int,Int}}}())
        init_constraints_idx_map(model,instance.id_con_idx_map)

        instance.get_num_scen = function()
            return num_scenarios(instance.internalModel)
        end
        instance.get_sense = function()
            return getObjectiveSense(instance.internalModel)
        end
        instance.get_status = function()
            return instance.status
        end
        instance.get_num_rows = function(id::Integer)
            return get_numcons(instance.internalModel,id)
        end
        instance.get_num_cols = function(id::Integer)
            return get_numvars(instance.internalModel,id)
        end
        instance.get_num_eq_cons = function(id::Integer)
            return length(instance.id_con_idx_map[id][1])
        end
        instance.get_num_ineq_cons = function(id::Integer)
            return length(instance.id_con_idx_map[id][2])
        end
        instance.set_status = function(s::Integer)
            instance.status = s
        end
        instance.set_num_rows = function(id::Integer,v::Integer) end
        instance.set_num_cols = function(id::Integer,v::Integer) end
        instance.set_num_eq_cons = function(id::Integer,v::Integer) end
        instance.set_num_ineq_cons = function(id::Integer,v::Integer) end


        instance.str_init_x0 = function(id,x0)
            assert(id in getScenarioIds(instance.internalModel))
            mm = get_model(instance.internalModel,id)
            nvar = get_numvars(instance.internalModel,id)
            @assert length(x0) == nvar

            for i=1:nvar
                x0[i] = getValue(Variable(mm,i))
                isnan(x0[i])?x0[i]=1.0:nothing
            end
        end  

        instance.str_prob_info = function(id,mode,clb,cub,rlb,rub)
            # @show id
            if mode == :Structure
                nn = get_numvars(instance.internalModel,id)
                mm = get_numcons(instance.internalModel,id)
                # @show nn,mm
                return (nn,mm)
            elseif mode == :Values
                # @show length(clb),length(cub)
                mm = get_model(instance.internalModel,id)
                nvar = get_numvars(instance.internalModel,id)
                @assert length(clb) == nvar 
                @assert length(cub) == nvar
                array_copy(mm.colUpper, 1, cub, 1, nvar)
                array_copy(mm.colLower, 1, clb, 1, nvar)
                lb,ub = getConstraintBounds(mm)
                (eq_idx, ieq_idx) = instance.id_con_idx_map[id]
                @assert length(lb) == length(ub)
                @assert length(eq_idx) + length(ieq_idx) == length(lb)
                num_eq = length(eq_idx)
                for i in eq_idx
                    rlb[i[2]] = lb[i[1]]
                    rub[i[2]] = ub[i[1]]
                end
                for i in ieq_idx
                    rlb[i[2]+num_eq] = lb[i[1]]
                    rub[i[2]+num_eq] = ub[i[1]]
                end
            else
                @assert false mode
            end
        end 


        instance.str_eval_f = function(id,x0,x1)
            e =  get_nlp_evaluator(instance.internalModel,id)
            return MathProgBase.eval_f(e,build_x(instance.internalModel,id,x0,x1))
        end

        instance.str_eval_g = function(id,x0,x1, new_eq_g, new_inq_g)
            e = get_nlp_evaluator(instance.internalModel,id)
            g = Vector{Float64}(get_numcons(instance.internalModel,id))
            MathProgBase.eval_g(e,g,build_x(instance.internalModel,id,x0,x1))
            (eq_idx, ieq_idx) = instance.id_con_idx_map[id]
            @assert length(new_eq_g) == length(eq_idx)
            @assert length(ieq_idx) == length(new_inq_g)
            for i in eq_idx
                new_eq_g[i[2]] = g[i[1]]
            end
            for i in ieq_idx
                new_inq_g[i[2]] = g[i[1]]
            end
        end

        instance.str_eval_grad_f = function(rowid, colid, x0, x1, new_grad_f)
            @assert rowid >= colid
            @assert sum(new_grad_f) == 0.0 
            e = get_nlp_evaluator(instance.internalModel,rowid)
            x = build_x(instance.internalModel,rowid,x0,x1)
            g = Vector{Float64}(length(x))
            MathProgBase.eval_grad_f(e,g,x)
            @assert length(g) == MathProgBase.numvar(get_model(instance.internalModel,rowid))
            @assert length(new_grad_f) == get_numvars(instance.internalModel,colid)

            var_idx_map = get_jac_col_var_idx(instance.internalModel,rowid,colid)
            for i in var_idx_map
                new_grad_f[i[2]] = g[i[1]]
            end
        end



        instance.str_eval_jac_g = function(rowid, colid, x0 , x1, mode, e_rowidx, e_colptr, e_values, i_rowidx, i_colptr, i_values)
            # @show "str_eval_jac_g", rowid, colid, mode
            m = instance.internalModel
            @assert rowid<=num_scenarios(m) && colid <= num_scenarios(m)
            @assert rowid >= colid
            if(mode == :Structure)
                e = get_nlp_evaluator(m,rowid)
                (jac_I,jac_J) = MathProgBase.jac_structure(e)
                # @show jac_I
                # @show jac_J
                (eq_idx, ieq_idx) = instance.id_con_idx_map[rowid]
                var_idx = get_jac_col_var_idx(m,rowid,colid)
                # @show var_idx

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

                # @show eq_jac_I
                # @show eq_jac_J
                eq_jac = sparse(eq_jac_I,eq_jac_J,ones(Float64,length(eq_jac_I)),length(eq_idx),get_numvars(m,colid))
                # @show ieq_jac_I
                # @show ieq_jac_J
                ieq_jac = sparse(ieq_jac_I,ieq_jac_J,ones(Float64,length(ieq_jac_I)),length(ieq_idx),get_numvars(m,colid))
                # @show (length(eq_jac.nzval), length(ieq_jac.nzval))

                return (length(eq_jac.nzval), length(ieq_jac.nzval))
            elseif(mode == :Values)
                e = get_nlp_evaluator(m,rowid)
                (jac_I,jac_J) = MathProgBase.jac_structure(e)
                jac_g = Vector{Float64}(length(jac_I))
                MathProgBase.eval_jac_g(e,jac_g,build_x(m,rowid,x0,x1))
                (eq_idx, ieq_idx) = instance.id_con_idx_map[rowid]
                var_idx = get_jac_col_var_idx(m,rowid,colid)

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
                    eq_jac = sparse(eq_jac_I,eq_jac_J,eq_jac_g, length(eq_idx),get_numvars(m,colid), keepzeros=true)
                    # @show eq_jac
                    # @show eq_jac.rowval
                    # @show eq_jac.colptr
                    # @show eq_jac.nzval
                    array_copy(eq_jac.rowval,1,e_rowidx,1,length(eq_jac.rowval))
                    array_copy(eq_jac.colptr,1,e_colptr,1,length(eq_jac.colptr))
                    array_copy(eq_jac.nzval, 1,e_values,1,length(eq_jac.nzval))
                    convert_to_c_idx(e_rowidx)
                    convert_to_c_idx(e_colptr)

                    filename = string("jaceq_",rowid,"_",colid)
                    write_mat_to_file(filename,eq_jac)
                end

                if(length(ieq_jac_g) != 0)
                    ieq_jac = sparse(ieq_jac_I,ieq_jac_J,ieq_jac_g, length(ieq_idx),get_numvars(m,colid), keepzeros=true)
                    # @show ieq_jac
                    # @show ieq_jac.rowval
                    # @show ieq_jac.colptr
                    # @show ieq_jac.nzval
                    array_copy(ieq_jac.rowval,1,i_rowidx,1,length(ieq_jac.rowval))
                    array_copy(ieq_jac.colptr,1,i_colptr,1,length(ieq_jac.colptr))
                    array_copy(ieq_jac.nzval, 1,i_values,1,length(ieq_jac.nzval))
                    convert_to_c_idx(i_rowidx)
                    convert_to_c_idx(i_colptr)
                    
                    filename = string("jacieq_",rowid,"_",colid)
                    write_mat_to_file(filename,ieq_jac)
                end
            else
                @assert false mode
            end
            # @show "end str_eval_jac_g"
        end


        instance.str_eval_h = function(rowid, colid, x0, x1, obj_factor, lambda, mode, rowidx, colptr, values)
            # @show "str_eval_h", rowid, colid, mode, length(x0), length(x1),obj_factor, length(lambda), length(rowidx), length(colptr), length(values)
            m = instance.internalModel
            @assert rowid<=num_scenarios(m) && colid <=num_scenarios(m)
            low = min(rowid, colid)
            high = max(rowid, colid)
            if(mode == :Structure)
                e = get_nlp_evaluator(m,high)
                (h_J, h_I) = MathProgBase.hesslag_structure(e) #reverse I, J as we want upper trangular matrix
                col_var_idx, row_var_idx = get_h_var_idx(m,rowid, colid)
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
                e = get_nlp_evaluator(m,high)
                (h_J, h_I) = MathProgBase.hesslag_structure(e)
                h = Vector{Float64}(length(h_I))
                # @show h_I
                # @show h_J
                # @show h
                MathProgBase.eval_hesslag(e,h,build_x(m,high,x0,x1),obj_factor,lambda)
                col_var_idx,row_var_idx = get_h_var_idx(m,rowid, colid)
                # @show col_var_idx
                # @show row_var_idx
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
                    (h0_I,h0_J) = MathProgBase.hesslag_structure(get_nlp_evaluator(m,0))
                    str_laghess = sparse([new_h_I;h0_I], [new_h_J;h0_J], [new_h;zeros(Float64,length(h0_I))], get_numvars(m,0),get_numvars(m,0), keepzeros=true)
                    # @show str_laghess.m, str_laghess.n, length(str_laghess.nzval)
                    # @show str_laghess
                    
                    array_copy(str_laghess.rowval,1,rowidx,1,length(str_laghess.rowval))
                    array_copy(str_laghess.colptr,1,colptr,1,length(str_laghess.colptr))
                    array_copy(str_laghess.nzval, 1,values,1,length(str_laghess.nzval))

                    # @show str_laghess.nzval
                else
                    str_laghess = sparse(new_h_I, new_h_J, new_h, get_numvars(m,colid), get_numvars(m,rowid), keepzeros=true)
                    # @show str_laghess.m, str_laghess.n, length(str_laghess.nzval)
                    # @show str_laghess
                    
                    array_copy(str_laghess.rowval,1,rowidx,1,length(str_laghess.rowval))
                    array_copy(str_laghess.colptr,1,colptr,1,length(str_laghess.colptr))
                    array_copy(str_laghess.nzval, 1,values,1,length(str_laghess.nzval))
                    # @show str_laghess.nzval
                end
                
                filename = string("hess_",rowid,"_",colid)
                write_mat_to_file(filename,str_laghess)

                convert_to_c_idx(rowidx)
                convert_to_c_idx(colptr)
            else
                @assert false mode
            end 
            # @show "end str_eval_h"
        end
        return instance
    end
end

#########################
# Helper 
#########################
function get_model(m,id)
    return id==0?m:getchildren(m)[id]
end

function get_numvars(m,id)
    mm = get_model(m,id)
    nvar = MathProgBase.numvar(mm) - length(getStructure(mm).othermap)
    return nvar
end

function get_numcons(m,id)
    mm = get_model(m,id)
    return MathProgBase.numconstr(mm)
end

function init_constraints_idx_map(m,map)
    assert(length(map) == 0)
    for id in getScenarioIds(m)
        eq_idx = Dict{Int,Int}()
        ieq_idx = Dict{Int,Int}()
        push!(map,id=>Pair(eq_idx,ieq_idx))
        lb,ub=getConstraintBounds(get_model(m,id))
        for i =1:length(lb)
            if lb[i] == ub[i]
                eq_idx[i] = length(eq_idx) + 1
            else
                ieq_idx[i] = length(ieq_idx) + 1 #remember to offset length(eq_idx)
            end
        end
    end
end

function get_nlp_evaluator(m,id)
    # @show id,getScenarioIds(m)
    # @show getProcIdxSet(m)
    # assert(id == 0 || id in getProcIdxSet(m))
    e = JuMPNLPEvaluator(get_model(m,id))
    MathProgBase.initialize(e,[:Grad,:Jac,:Hess])
    return e
end

function build_x(m,id,x0,x1)
    # @show "build_x", id, length(x0), length(x1)
    # @show x0, x1
    if id==0
        return x0
    else
        #build x using index tracking info
        mm = get_model(m,id)
        othermap = getStructure(mm).othermap
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


function get_jac_col_var_idx(m,rowid, colid)  #this method customerized for no linking constraint presented
    # @show "get_jac_col_var_idx",rowid,colid
    idx_map = Dict{Int,Int}() #dummy (jump) -> actual used
    if rowid == colid
        nvar = get_numvars(m,rowid)
        for i = 1:nvar
            idx_map[i] = i
        end
    else
        @assert rowid!=0 && rowid != colid
        mm = get_model(m,rowid)
        othermap = getStructure(mm).othermap
        for p in othermap
            pidx = p[1].col
            cidx = p[2].col
            idx_map[cidx] = pidx
        end
    end
    # @show idx_map
    return idx_map
end

function get_h_var_idx(m,rowid, colid)
    # @show "get_h_var_idx",rowid,colid
    col_idx_map = Dict{Int,Int}() #dummy (jump) -> actual used
    row_idx_map = Dict{Int,Int}()
    if rowid == colid
        nvar = get_numvars(m,rowid)
        for i = 1:nvar
            col_idx_map[i] = i
            row_idx_map[i] = i
        end
    elseif rowid == 0  #border
        mm = get_model(m,colid)
        othermap = getStructure(mm).othermap
        for p in othermap
            pidx = p[1].col
            cidx = p[2].col
            col_idx_map[cidx] = pidx
        end
        for i = 1:get_numvars(m,colid)
            row_idx_map[i] = i
        end
    elseif colid == 0 #root contrib.
        mm = get_model(m,rowid)
        othermap = getStructure(mm).othermap
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

        for c = 1:N
            crange = nzrange(actual,c)
            full.nzval[crange] = actual.nzval[crange] 
        end        
        return full
    end
end


# function array_copy(src,dest)
#     @assert length(src)==length(dest)
#     array_copy(src,1,dest,1,length(src))
# end

function array_copy(src,os, dest, od, n)
    for i=0:n-1
        dest[i+od] = src[i+os]
    end
end

function write_mat_to_file(filename,mat)
    if(false)
        filename = string("./mat/",filename)
        pre_filename = string(filename,"_0")
        i = 0
        while isfile(pre_filename)
            i += 1
            pre_filename = string(filename,"_",i)
        end
        @show "output : ", pre_filename
        writedlm(pre_filename,mat,",")
    end
end

function convert_to_c_idx(indicies)
    for i in 1:length(indicies)
        indicies[i] = indicies[i] - 1
    end
end

#######
# Linking with PIPS Julia Structure interface
######
function createStructJuMPPipsNlpProblem(model::JuMP.Model)
    @show "createStructJuMPPipsNlpProblem"
    comm = MPI.COMM_WORLD
    @show "[$(MPI.Comm_rank(comm))/$(MPI.Comm_size(comm))] create problem "
    @show comm

    prob = createProblemStruct(comm, StructJuMPModel(model,0))
    @show "end createStructJuMPPipsNlpProblem"
    return prob
end

function solveStructJuMPPipsNlpProblem(prob)
    @show "solveStructJuMPPipsNlpProblem"
    @show prob
    ret = solveProblemStruct(prob)
    @show "end solveStructJuMPPipsNlpProblem"
    return ret
end

function solve(model)
    @show "solve"
    global m = model
    MPI.Init()

    prob = createStructJuMPPipsNlpProblem(model)

    solveStructJuMPPipsNlpProblem(prob)

    freeProblemStruct(prob)

    MPI.Finalize()
end

end

