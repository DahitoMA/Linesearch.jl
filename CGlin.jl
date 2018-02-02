# A version of the conjugate gradient method for linesearch.

"""A version of the conjugate gradient method for linesearch.
CG(A, b, ϵa, ϵr, itmax; quad) solves the symmetric linear system 'A * x = b'
If quad = true, the values of the quadratic model are computed
A can be positive definite or not.
"""
function CGlin(A, b, ϵa::Float64=1e-8, ϵr::Float64=1e-6, itmax::Int=0; quad::Bool=false)
    n = size(b, 1) # size of the problem
    (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size")
    @info(loggerCGlin, @sprintf("CGlin: system of %d equations in %d variables", n, n))

    x = zeros(n) # initial estimation x = 0
    xNorm = 0.0
    xNorms = [xNorm] # Values of ‖x‖
    r = -b # initial residual r = Ax-b = -b
    d = b # first descent direction
    rNorm = norm(r, 2)
    ϵ = ϵa + ϵr * rNorm
    if quad
        q = 0.0
        qvalues = [q] # values of the quadratic model
        @info(loggerCGlin, @sprintf("%5s %10s %10s\n", "Iter", "‖r‖", "q"))
        @info(loggerCGlin, @sprintf("    %d    %8.1e    %8.1e", iter, rNorm, q))
    end

    iter = 0
    itmax == 0 && (itmax = 2 * n)

    solved = rNorm ≤ ϵ
    tired = iter ≥ itmax

    while ! (solved || tired)
        iter += 1
        Ad = A * d
        dAd = dot(d, Ad)

        # if the model is not convexe, the algorithm stops
        if dAd ≤ ϵ * dot(d, d)
            @debug(loggerCGlin, @sprintf("non positive curvature dAd = %8.1e", dAd))
            iter == 1 && return b
            return x
        end

        α = rNorm^2 / dAd # step for x estimation
        x = x + α * d # new estimation
        xNorm = norm(x, 2)
        push!(xNorms, xNorm)
        roldNorm = rNorm
        r = r + α * Ad # new residual
        rNorm = norm(r, 2)

        if quad
            q = -dot(b, x) + 0.5 * dot(x, A * x)
            push!(qvalues, q)
            @info(loggerCGlin, @sprintf("    %d    %8.1e    %8.1e", iter, rNorm, q))
        end

        solved = rNorm ≤ ϵ
        tired = iter ≥ itmax

        (solved || tired) && continue
        β = rNorm^2 / roldNorm^2 # step for the next descent direction
        d = -r + β * d # new descent direction

        end
        return x
    end
