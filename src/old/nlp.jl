# Add parent::Int here to avoid overwriting same method in JuMP. Should be safe.
function JuMP.parseNLExpr_runtime(m::JuMP.Model, x::JuMP.Variable, tape, parent::Int, values)
    #JuMP.__last_model[1] = x.m  #__last_model is not defined in JuMP
    if x.m === m
        push!(tape, JuMP.NodeData(ReverseDiffSparse.VARIABLE, x.col, parent))
    else
        othermap = getStructure(m).othermap
        if haskey(othermap, x)
            newx = othermap[x]
        else
            newx = JuMP.Variable(m,JuMP.getlowerbound(x),JuMP.getupperbound(x),:Cont, JuMP.getname(x))
            othermap[x] = newx  
        end
        push!(tape, JuMP.NodeData(ReverseDiffSparse.VARIABLE, newx.col, parent))
    end
    nothing
end