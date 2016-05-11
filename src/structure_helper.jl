# structure_helper.jl
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

function get_var_value(m,id)
    mm = get_model(m,id)
    v = Float64[];
    for i = 1:get_numvars(m,id)
        v = [v;getValue(Variable(mm,i))]
    end
    return v
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
                ieq_idx[i] =length(ieq_idx) + 1 #remember to offset length(eq_idx)
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
        # @show x0
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
        # @show new_x
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
        nvar = get_numvars(m,rowid) #need to place model variable in front of non model variable.
        for i = 1:nvar
            col_idx_map[i] = i
            row_idx_map[i] = i
        end
    elseif rowid == 0  && colid != 0 #border
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
    elseif colid == 0 && rowid != 0 #root contrib.
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
            for i=nzrange(actual,c)
                r = actual.rowval[i]
                v = actual.nzval[i]
                if(v!=0)
                    full[r,c] = v
                end
            end  
            # full.nzval[crange] = actual.nzval[crange] 
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

### for non structure solver

function strip_x(m,id,x,start_idx)
    mm = get_model(m,id)
    nx = get_numvars(m,id)
    new_x = Vector{Float64}(MathProgBase.numvar(mm))
    array_copy(x,start_idx,new_x,1,nx)

    othermap = getStructure(mm).othermap
    for i in othermap
        pid = i[1].col
        cid = i[2].col
        new_x[cid] = x[pid]
        @assert cid > nx
    end
    return new_x
end

function g_numcons(m)
    ncon = 0
    for i=0:num_scenarios(m)
        ncon += get_numcons(m,i)
    end
    return ncon
end

function g_numvars(m)
    nvar = 0
    for i=0:num_scenarios(m)
        nvar += get_numvars(m,i)
    end
    return nvar
end
