# Implementation of Armijo backtracking linesearch
# Armijo(model, x, s, fx, g, t, α) computes a step t satisfying Armijo condition
# and returns the new iterate xtrial

"""Armijo backtracking linesearch
The function returns the new iterate for the linesearch
"""
function Armijo(model, x, s, fx=obj(model, x), g=grad(model, x), t=1., α=0.8)
    @info(loggerArm, @sprintf("Armijo backtracking linesearch on problem %s", string(model)))
    (t ≤ 0) && (t = 1.)
    @debug(loggerArm, @sprintf("initial t = %7.1e", t))
    (α ≤ 0 || α ≥ 1) && (α = 0.8)
    @debug(loggerArm, @sprintf("α = %8.1e", α))
    xtrial = x + t * s
    fxtrial = obj(model, xtrial)
    αgs = α * dot(g, s)
    β = fx + t * αgs

    while fxtrial > β && t > eps()
        @debug(loggerArm, @sprintf("fxtrial = %8.1e > β = %8.1e", fxtrial, β))
        t = 0.5 * t
        @debug(loggerArm, @sprintf("new t = %7.1e", t))
        xtrial = x + t * s
        fxtrial = obj(model, xtrial)
        β = fx + t * αgs
    end
    return xtrial
end
