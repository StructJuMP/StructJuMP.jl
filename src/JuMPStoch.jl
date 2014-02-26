module JuMPStoch

# import JuMP.JuMPDict
importall JuMP

using MathProgBase
using MathProgBase.MathProgSolverInterface

importall Base

using Base.Meta

export StochasticData, StochasticModel, getStochastic, StochasticBlock, ancestor, @defStochasticVar

# JuMP rexports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint,
# Functions
    # Relevant to all
    print,show,
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addVar, addVars, solve, copy,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    # Macros and support functions
    @addConstraint, @defVar, 
    @defConstrRef, @setObjective, addToExpression

pushchild!(m::Model, block) = push!(m.ext[:Stochastic].children, block)

type StochasticData
    id
    children::Vector{Model}
    parent
    vars::Dict{String,Variable}
end

StochasticData() = StochasticData(nothing,Model[],nothing,Dict{String,Variable}())

function StochasticModel(;solver=nothing)
    m = Model(solver=solver)
    m.ext[:Stochastic] = StochasticData()
    return m
end

function StochasticModel(id, children, parent)
    m = Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(id, children, parent,Dict{String,Variable}())
    return m
end

function getStochastic(m::Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

function StochasticBlock(m::Model, id)
    stoch = getStochastic(m)
    ch = StochasticModel(id, Model[], m)
    pushchild!(m, ch)
    return ch
end

function Variable(m::Model,lower::Number,upper::Number,cat::Int,name::String)
    m.numCols += 1
    push!(m.colNames, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,NaN)
    var = Variable(m, m.numCols)
    stoch = getStochastic(m)
    stoch.vars[string(name)] = var
    return var
end


function ancestor(m::Model, level::Int)
    stoch = getStochastic(m)
    if level == 1
        return stoch.vars
    elseif level > 1
        return ancestor(stoch.parent, level-1)
    else
        error("Can only treat positive levels")
    end
end
ancestor(m::Model) = ancestor(m, 1)

macro defStochasticVar(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = -Inf
        ub = Inf
    end
    t = JuMP.CONTINUOUS
    if length(extra) > 0
        gottype = 0
        if extra[1] == :Int || extra[1] == :Bin
            gottype = 1
            if extra[1] == :Int
                t = JuMP.INTEGER
            else
                if lb != -Inf || ub != Inf
                    error("Bounds may not be specified for binary variables. These are always taken to have a lower bound of 0 and upper bound of 1.")
                end
                t = JuMP.INTEGER
                lb = 0.0
                ub = 1.0
            end
        end
        if length(extra) - gottype == 3
            # adding variable to existing constraints
            objcoef = esc(extra[1+gottype])
            cols = esc(extra[2+gottype])
            coeffs = esc(extra[3+gottype])
            if !isa(var,Symbol)
                error("Cannot create multiple variables when adding to existing constraints")
            end
            return quote
                $(esc(var)) = Variable($m,$lb,$ub,$t,$objcoef,$cols,$coeffs,name=$(string(var)))
                nothing
            end
        elseif length(extra) - gottype != 0
            error("Syntax error in defVar")
        end
    end

    #println("lb: $lb ub: $ub var: $var")      
    if isa(var,Symbol)
        # easy case
        return quote
            $(esc(var)) = Variable($m,$lb,$ub,$t,$(string(var)))
            nothing
        end
    else
        if !isexpr(var,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end
        varname = esc(var.args[1])
        idxvars = {}
        idxsets = {}
        refcall = Expr(:ref,varname)
        for s in var.args[2:end]
            if isa(s,Expr) && s.head == :(=)
                idxvar = s.args[1]
                idxset = esc(s.args[2])
            else
                idxvar = gensym()
                idxset = esc(s)
            end
            push!(idxvars, idxvar)
            push!(idxsets, idxset)
            push!(refcall.args, esc(idxvar))
        end
        tup = Expr(:tuple, [esc(x) for x in idxvars]...)
        code = :( $(refcall) = Variable($m, $lb, $ub, $t, $(string(var.args[1]))*string($tup) ) )
        # code = :( $(refcall) = Variable($m, $lb, $ub, $t) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $(esc(idxvar)) in $idxset
                    $code
                end
            end
        end
       
        mac = Expr(:macrocall,symbol("@gendict"),varname,:Variable,idxsets...)
        addDict = :( push!($(m).dictList, $varname) )
        code = quote 
            $mac
            $code
            $addDict
            nothing
        end
        return code
    end
end

end
