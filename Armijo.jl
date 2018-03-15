# Implementation of Armijo backtracking linesearch

"""Armijo backtracking linesearch
Armijo(model, x, s, fx, g, t, α) computes a step t satisfying Armijo condition
The function returns (xtrial, fxtrial) : the new iterate for the linesearch and
the value of the objective in xtrial
"""
function Armijo(model, x, s, fx=obj(model, x), g=grad(model, x), t=1., α=1e-4)
    @info(loggerArm, @sprintf("Armijo backtracking linesearch on problem %s", model.meta.name))
    (t ≤ 0) && (t = 1.)
    @debug(loggerArm, @sprintf("initial step t = %7.1e", t))
    (α ≤ 0 || α ≥ 1) && (α = 1e-4)
    @debug(loggerArm, @sprintf("α = %8.1e", α))
    xtrial = x + t * s
    fxtrial = obj(model, xtrial) # f(x + t * s)
    @debug(loggerArm, @sprintf("f(x) = %8.1e", fxtrial))
    αgs = α * dot(g, s)
    β = fx + t * αgs # f(x) + t * α * gᵀs
    @debug(loggerArm, @sprintf("β = %8.1e", β))

    while fxtrial > β && t > eps()
        @debug(loggerArm, @sprintf("f(x) = %8.1e > β = %8.1e", fxtrial, β))
        t = 0.5 * t
        @debug(loggerArm, @sprintf("new t = %7.1e", t))
        xtrial = x + t * s
        fxtrial = obj(model, xtrial)
        β = fx + t * αgs
    end
    @debug(loggerArm, @sprintf("final values : f(x+t*s) = %8.1e ≤ β = %8.1e", fxtrial, β))
    return xtrial, fxtrial
end
