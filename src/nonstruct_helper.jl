# nonstruct_helper.jl

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