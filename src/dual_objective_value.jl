const MOIU = MOI.Utilities
using LinearAlgebra

# Taken from https://github.com/JuliaOpt/MathOptInterface.jl/pull/473
"""
    set_dot(x::Vector, y::Vector, set::AbstractVectorSet)

Return the scalar product between a vector `x` of the set `set` and a vector
`y` of the dual of the set `s`.
"""
set_dot(x::Vector, y::Vector, set::MOI.AbstractVectorSet) = dot(x, y)

"""
    set_dot(x, y, set::AbstractScalarSet)

Return the scalar product between a number `x` of the set `set` and a number
`y` of the dual of the set `s`.
"""
set_dot(x, y, set::MOI.AbstractScalarSet) = dot(x, y)

_getconstant(set::MOI.AbstractScalarSet, T) = MOIU.getconstant(set)
_getconstant(set::MOI.Integer, T) = zero(T)


scalar_constant(T::Type, ::MOI.SingleVariable) = zero(T)
scalar_constant(::Type, f::MOI.AbstractScalarFunction) = MOI._constant(f)

function constraint_constant(model::MOI.ModelLike,
                             ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,
                                                     <:MOI.AbstractVectorSet},
                             T::Type)
    return MOIU._constant(MOI.get(model, MOI.ConstraintFunction(), ci))
end
function constraint_constant(model::MOI.ModelLike,
                             ci::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction,
                                                     <:MOI.AbstractScalarSet},
                             T::Type)
    return scalar_constant(T, MOI.get(model, MOI.ConstraintFunction(), ci)) -
           _getconstant(MOI.get(model, MOI.ConstraintSet(), ci), T)
end

"""
    objective_bound(model::MOI.ModelLike,
                    F::Type{<:MOI.AbstractFunction},
                    S::Type{<:MOI.AbstractSet},
                    T::Type)

Return the part of `ObjectiveBound` due to the constraint of index `ci` using
scalar type `T`.
"""
function objective_bound(model::MOI.ModelLike,
                         ci::MOI.ConstraintIndex,
                         T::Type)
    return set_dot(constraint_constant(model, ci, T),
                   MOI.get(model, MOI.ConstraintDual(), ci),
                   MOI.get(model, MOI.ConstraintSet(), ci))
end

function objective_bound(model::MOI.ModelLike,
                         ci::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction,
                                                 <:MOI.Interval},
                         T::Type)
    constant = scalar_constant(T, MOI.get(model, MOI.ConstraintFunction(), ci))
    set = MOI.get(model, MOI.ConstraintSet(), ci)
    dual = MOI.get(model, MOI.ConstraintDual(), ci)
    if dual < zero(dual)
        # The dual is negative so it is in the dual of the MOI.LessThan cone
        # hence the upper bound of the Interval set is tight
        constant -= set.upper
    else
        # the lower bound is tight
        constant -= set.lower
    end
    return set_dot(constant, dual, set)
end

"""
    objective_bound(model::MOI.ModelLike,
                    F::Type{<:MOI.AbstractFunction},
                    S::Type{<:MOI.AbstractSet},
                    T::Type)

Return the part of `ObjectiveBound` due to `F`-in-`S` constraints using scalar
type `T`.
"""
function objective_bound(model::MOI.ModelLike,
                         F::Type{<:MOI.AbstractFunction},
                         S::Type{<:MOI.AbstractSet},
                         T::Type)
    bound = zero(T) # sum won't work if there are now constraints
    for ci in MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
        bound += objective_bound(model, ci, T)
    end
    return bound
end

function objective_bound(model::MOI.ModelLike,
                         F::Type{MOI.VectorOfVariables},
                         S::Type{<:MOI.AbstractVectorSet},
                         T::Type)
    # No constant in the function nor set so no contribution to the objective
    # bound
    return zero(T)
end


"""
    dual_objective_value(model::MOI.ModelLike, ::MOI.ObjectiveBound, T::Type)::T

Compute the objective value of the dual problem of type `T` using the
`ConstraintDual` results and the `ConstraintFunction` and `ConstraintSet` values.
"""
function dual_objective_value(model::MOI.ModelLike, T::Type)
    bound = zero(T) # sum will not work if there are zero constraints
    for (F, S) in MOI.get(model, MOI.ListOfConstraints())
        bound += objective_bound(model, F, S, T)::T
    end
    if MOI.get(model, MOI.ObjectiveSense()) != MOI.MAX_SENSE
        bound = -bound
    end
    dual_status = MOI.get(model, MOI.DualStatus())
    if dual_status == MOI.INFEASIBILITY_CERTIFICATE ||
        status == MOI.NEARLY_INFEASIBILITY_CERTIFICATE
        # The objective constant should not be present in rays
        F = MOI.get(model, MOI.ObjectiveFunctionType())
        f = MOI.get(model, MOI.ObjectiveFunction{F}())
        bound += scalar_constant(T, f)
    end
    return bound::T
end
