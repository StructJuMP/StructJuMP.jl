# Add parent::Int here to avoid overwritting same method in JuMP. Should be safe.
function JuMP.parseNLExpr_runtime(m::JuMP.Model, x::JuMP.Variable, tape, parent::Int, values)
    # @show m
    # @show parent
    # @show JuMP.getName(x)
    # @show x.col
    # @show tape
    # @show parent
    # @show values

    JuMP.__last_model[1] = x.m
    if x.m === m
        push!(tape, JuMP.NodeData(ReverseDiffSparse.VARIABLE, x.col, parent))
    else
        # @show "othervars"
        # others = getStochastic(m).othervars
        othermap = getStochastic(m).othermap
        # push!(others, x)
        if haskey(othermap, x)
            newx = othermap[x]
        else
            newx = JuMP.Variable(m,JuMP.getLower(x),JuMP.getUpper(x),:Cont, JuMP.getName(x))
            othermap[x] = newx
        end
        push!(tape, JuMP.NodeData(ReverseDiffSparse.VARIABLE, newx.col, parent))
    end
    nothing
end

