module JuMPStoch

import JuMP.JuMPDict
using JuMP

using MathProgBase
using MathProgBase.MathProgSolverInterface

importall Base

using Base.Meta

export StochasticData, StochasticModel, getStochastic, parent, children, StochasticBlock, StochasticVariable, variables, @defStochasticVar

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

pushchild!(m::Model, block) = push!(getStochastic(m).children, block)

type StochasticData
    id
    children::Vector{Model}
    parent
    varstup::Dict{Tuple,Variable}
    vardict::Dict
end

StochasticData() = StochasticData(nothing,Model[],nothing,Dict{Tuple,Variable}(),Dict())

function StochasticModel(;solver=nothing)
    m = Model(solver=solver)
    m.ext[:Stochastic] = StochasticData()
    return m
end

function StochasticModel(id, children, parent)
    m = Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(id, children, parent,Dict{Tuple,Variable}(),Dict())
    return m
end

function getStochastic(m::Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

parent(m::Model) = getStochastic(m).parent
children(m::Model) = getStochastic(m).children

function StochasticBlock(m::Model, id)
    stoch = getStochastic(m)
    ch = StochasticModel(id, Model[], m)
    pushchild!(m, ch)
    return ch
end

function StochasticVariable(m::Model,lower::Number,upper::Number,cat::Int,name::String)
    m.numCols += 1
    push!(m.colNames, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,NaN)
    var = Variable(m, m.numCols)
    stoch = getStochastic(m)
    stoch.varstup[tuple(name)] = var
    return var
end

function StochasticVariable(m::Model,lower::Number,upper::Number,cat::Int,args::Tuple)
    m.numCols += 1
    push!(m.colNames, "")
    # push!(m.colNames, name)
    push!(m.colLower, convert(Float64,lower))
    push!(m.colUpper, convert(Float64,upper))
    push!(m.colCat, cat)
    push!(m.colVal,NaN)
    var = Variable(m, m.numCols)
    stoch = getStochastic(m)
    stoch.varstup[args] = var
    return var
end

variables(m::Model) = getStochastic(m).vardict

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
                $(esc(var)) = StochasticVariable($m,$lb,$ub,$t,$objcoef,$cols,$coeffs,name=$(string(var)))
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
            $(esc(var)) = StochasticVariable($m,$lb,$ub,$t,$(string(var)))
            $(m).ext[:Stochastic].vardict[$(quot(var))] = $(esc(var))   
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
        code = :( $(refcall) = StochasticVariable($m, $lb, $ub, $t, tuple($(string(var.args[1])), $(tup)...)) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $(esc(idxvar)) in $idxset
                    $code
                end
            end
        end

        mac = Expr(:macrocall,symbol("@genStochDict"),varname,:Variable,idxsets...)
        addVarDict = :( $(m).ext[:Stochastic].vardict[$(quot(var.args[1]))] = $varname )   
        addDict = :( push!($(m).dictList, $varname) )
        code = quote
            $mac
            $code
            $addDict
            $addVarDict
            nothing
        end
        return code
    end
end

getindex(d::Variable, owner::Model) = getStochastic(owner).varstup[tuple(string(d))]
getindex(d::JuMPDict,owner::Model,args...) = getStochastic(owner).varstup[tuple(d.name,args...)]

macro genStochDict(instancename,T,idxsets...)
    N = length(idxsets)
    typename = symbol(string("JuMPDict",gensym()))
    isrange = Array(Bool,N)
    offset = Array(Int,N)
    dictnames = Array(Symbol,N)
    for i in 1:N
        if isexpr(idxsets[i],:(:)) && length(idxsets[i].args) == 2 # don't yet optimize ranges with steps
            isrange[i] = true
            if isa(idxsets[i].args[1],Int)
                offset[i] = 1 - idxsets[i].args[1]
            else
                error("Currently only ranges with integer compile-time starting values are allowed as index sets. $(idxsets[i].args[1]) is not an integer in range $(idxsets[i]).")
            end
        else
            isrange[i] = false
            dictnames[i] = gensym()
        end
    end
    typecode = :(type $(typename){T} <: JuMPDict{T}; innerArray::Array{T,$N}; name::String;
                        indexsets end)
    builddicts = quote end
    for i in 1:N
        if !isrange[i]
            push!(typecode.args[3].args,:($(dictnames[i])::Dict))
            push!(builddicts.args, quote
                $(esc(dictnames[i])) = Dict();
                for (j,k) in enumerate($(esc(idxsets[i])))
                    $(esc(dictnames[i]))[k] = j
                end
            end)
        end
    end
    getidxlhs = :(getindex(d::$(typename)))
    setidxlhs = :(setindex!(d::$(typename),val))
    getidxrhs = :(getindex(d.innerArray))
    setidxrhs = :(setindex!(d.innerArray,val))
    getidxlhs2 = :(getindex(d::$(typename),owner::Model))
    getidxrhs2 = :(error("Must index JuMPDict"))
    maplhs = :(mapvals(f,d::$(typename)))
    maprhs = :($(typename)(map(f,d.innerArray),d.name,d.indexsets))
    for i in 1:N
        varname = symbol(string("x",i))

        if isrange[i]
            push!(getidxlhs.args,:($varname))
            push!(setidxlhs.args,:($varname))

            push!(getidxrhs.args,:($varname+$(offset[i])))
            push!(setidxrhs.args,:($varname+$(offset[i])))
        else
            push!(getidxlhs.args,varname)
            push!(setidxlhs.args,varname)

            push!(getidxrhs.args,:(d.($(Expr(:quote,dictnames[i])))[$varname]))
            push!(setidxrhs.args,:(d.($(Expr(:quote,dictnames[i])))[$varname]))
            push!(maprhs.args,:(d.($(Expr(:quote,dictnames[i])))))
        end
    end

    funcs = :($getidxlhs2 = $getidxrhs2; $getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs; $maplhs = $maprhs)
    geninstance = :($(esc(instancename)) = $(typename)(Array($T),$(string(instancename)),$(esc(Expr(:tuple,idxsets...)))))
    for i in 1:N
        push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
        if !isrange[i]
            push!(geninstance.args[2].args, esc(dictnames[i]))
        end
    end
    eval(Expr(:toplevel, typecode))
    eval(Expr(:toplevel, funcs))

    quote
        $builddicts
        $geninstance
    end

end

end
