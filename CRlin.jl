# A version of Stiefel’s Conjugate Residual method for linesearch.

"""A version of Stiefel’s Conjugate Residual method for linesearch.
CR(A, b, ϵa, ϵr, itmax; quad) solves the symmetric linear system 'A * x = b'
or the least-squares problem : 'min ‖b - A * x‖²'
If quad = true, the values of the quadratic model are computed
A can be positive definite or not.
"""
function CRlin(A, b, ϵa::Float64=1e-8, ϵr::Float64=1e-6, itmax::Int=0, ε::Float64=1e-6; quad::Bool=false)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(loggerCRlin, @sprintf("CRlin: system of %d equations in %d variables", n, n))

    x = zeros(n) # initial estimation x = 0
    xNorm = 0.0
    xNorms = [xNorm] # Values of ‖x‖
    r = b # initial residual r = b - Ax = b
    rNorm = norm(r, 2) # ‖r‖
    rNorm² = rNorm * rNorm
    pr = rNorm²
    s = A * r
    ρ = dot(r, s)
    p = r
    pNorm² = rNorm²
    q = s
    ϵ = ϵa + ϵr * rNorm
    pAp = ρ # = dot(p, q) = dot(r, s)

    iter = 0
    itmax == 0 && (itmax = 2 * n)

    if quad
        m = 0.0
        mvalues = [m] # values of the quadratic model
        @info(loggerCRlin, @sprintf("%5s %7s %8s", "Iter", "‖r‖", "q"))
        @info(loggerCRlin, @sprintf("%5d %7.1e %8.1e", iter, rNorm, m))
    end

    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax

    while ! (solved || tired)
        iter += 1

        if (pAp ≤ ε * pNorm²) || (ρ ≤ ε * rNorm²)
            @debug(loggerCRlin, @sprintf("nonpositive curvature detected: pAp = %8.1e and rAr = %8.1e", pAp, ρ))
            iter == 1 && return b
            return x
        end

        α = ρ / dot(q, q) # step
        x = x + α * p
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        r = r - α * q  # residual
        rNorm² = abs(rNorm² - α * ρ)
        rNorm = sqrt(rNorm²)

        if quad
            m = -dot(b, x) + 0.5 * dot(x, A * x)
            push!(mvalues, m)
            @info(loggerCRlin, @sprintf("%5d %7.1e %8.1e", iter, rNorm, m))
        end

        solved = rNorm <= ϵ
        tired = iter >= itmax
        (solved || tired) && continue

        s = A * r
        ρbar = ρ
        ρ = dot(r, s)
        β = ρ / ρbar # step for the direction calculus
        p = r + β * p # descent direction
        pNorm² = rNorm² + 2 * β * pr - 2 * β * α * pAp + β^2 * pNorm²
        pr = rNorm² + β * pr - β * α * pAp # pᵀr
        q = s + β * q
        pAp = ρ + β^2 * pAp # dot(p, q)

    end

    return x
end
