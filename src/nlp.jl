# Add parent::Int here to avoid overwritting same method in JuMP. Should be safe.
function JuMP.parseNLExpr_runtime(m::JuMP.Model, x::JuMP.Variable, tape, parent::Int, values)
    JuMP.__last_model[1] = x.m
    if x.m === m
        push!(tape, JuMP.NodeData(ReverseDiffSparse.VARIABLE, x.col, parent))
    else
        others = getStochastic(m).othervars
        push!(others, x)
        push!(tape, JuMP.NodeData(ReverseDiffSparse.EXTRA, length(others), parent))
    end
    nothing
end
