# structure_helper.jl

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


# function array_copy(src,dest)
#     @assert length(src)==length(dest)
#     array_copy(src,1,dest,1,length(src))
# end

