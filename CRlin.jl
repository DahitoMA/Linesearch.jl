import Krylov

# A version of Stiefel’s Conjugate Residual method for linesearch.

"""A version of Stiefel’s Conjugate Residual method for linesearch.
CR(A, b, ϵa, ϵr, itmax) solves the symmetric linear system 'A * x = b'
or the least-squares problem : 'min ‖b - A * x‖²'
A can be positive definite or not.
"""
function CRlin(A, b, ϵa::Float64=1e-8, ϵr::Float64=1e-6, itmax::Int=0; args...)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(loggerCRlin, @sprintf("CRlin: system of %d equations in %d variables", n, n))

    x = zeros(n) # initial estimation x = 0
    xNorm = 0.0
    xNorms = [xNorm] # Values of ‖x‖
    r = b # initial residual r = b - Ax = b
    rNorm = norm(r, 2) # ‖r‖
    rNorm² = rNorm * rNorm
    s = A * r
    ρ = dot(r, s)
    p = r
    q = s
    m = 0.0
    mvalues = [m] # values of the quadratic model
    ϵ = ϵa + ϵr * rNorm
    pAp = ρ # = dot(p, q) = dot(r, s)

    iter = 0
    itmax == 0 && (itmax = 2 * n)
    @info(loggerCRlin, @sprintf("%5s %7s %8s", "Iter", "‖r‖", "q"))
    @info(loggerCRlin, @sprintf("%5d %7.1e %8.1e", iter, rNorm, m))

    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax
    on_boundary = false

    while ! (solved || tired)
        iter += 1

        if pAp ≤ 0 || ρ ≤ 0
            @debug(loggerCRlin, @sprintf("nonpositive curvature detected: pAp = %8.1e and rAr = %8.1e", pAp, ρ))
            iter == 1 && return b
            return x
        end

        α = ρ / dot(q, q) # step
        x = x + α * p
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        Ax = A * x
        m = -dot(b, x) + 0.5 * dot(x, Ax)
        push!(mvalues, m)
        r = r - α * q  # residual
        rNorm² = abs(rNorm² - α * ρ)
        rNorm = sqrt(rNorm²)

        @info(loggerCRlin, @sprintf("%5d %7.1e %8.1e", iter, rNorm, m))

        solved = rNorm <= ϵ
        tired = iter >= itmax
        (solved || tired) && continue

        s = A * r
        ρbar = ρ
        ρ = dot(r, s)
        β = ρ / ρbar # step for the direction calculus
        p = r + β * p # descent direction
        q = s + β * q
        pAp = ρ + β^2 * pAp # dot(p, q)

    end

    return x
end
