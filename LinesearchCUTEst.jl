# An implementation of the linesearch method

"""An implementation of the linesearch method.
Linesearch(model, algo ; ϵa, ϵr, itemax) solves a continuous optimization
problem 'model' with absolute and relative tolerances ϵa and ϵr.
Steps are calculated using the argument 'algo', an optimization method using linesearch
"""
function LinesearchCUTEst(model, algo ; ϵa::Float64=1e-6, ϵr::Float64=1e-6, itemax::Int=10000)

    @info(loggerLinCUTEst, @sprintf("Linesearch: resolution of %s using %s", model.meta.name, string(algo)))
    x = model.meta.x0 # initial estimation from the model
    n = model.meta.nvar # size of the problem
    g = grad(model, x) # ∇f(x_0)
    H = hess_op(model, x) # ∇²f(x_0)

    normg0 = norm(g,2) # ‖∇f(x_0)‖
    normg = normg0 # ‖∇f(x_k)‖
    ϵ = ϵa + ϵr * normg # tolerance

    fx0 = obj(model, x) # f(x_0)
    fx = fx0 # f(x_k)

    fvalues = [fx] # improving values of the objective
    k = 0 # number of iterations

    @debug(loggerLinCUTEst, @sprintf("%4s  %9s  %7s", "Iter", "f", "‖∇f‖"))
    @debug(loggerLinCUTEst, @sprintf("%4d  %9.2e  %7.1e", k, fx, normg))

    while normg > ϵ && k < itemax # stopping criterion : ‖∇f(x_k)‖ <= ϵ or k >= itemax
        k += 1
        s = algo(H, -g, 0., min(0.1, sqrt(normg)), quad=true) # descent direction
        X = Armijo(model, x, s, fx, g)
        x = X[1]
        fx = X[2]
        fvalues = push!(fvalues, fx)
        g = grad(model, x)
        normg = norm(g,2)
        @debug(loggerLinCUTEst, @sprintf("%4d  %8.1e  %7.1e", k, fx, normg))
        H = hess_op(model, x)
    end

    @info(loggerLinCUTEst, @sprintf("%30s %s %9s %9s %9s %9s %3s %3s %4s %4s","name", "nvar", "f(x*)", "/ f(x0)", "‖∇f(x*)‖", "/ ‖∇f(x0)‖", "#f", "#g", "#Hv", "#it"))
    @info(loggerLinCUTEst, @sprintf("%30s %d  %8.1e  %8.1e    %7.1e  %7.1e     %d   %d   %d    %d", model.meta.name, n, fx, fx0, normg, normg0, neval_obj(model), neval_grad(model), neval_hprod(model), k))

    optimal = normg ≤ ϵ
    tired = k ≥ itemax
    status = optimal ? "optimal" : "tired"

    fx = @sprintf("%8.3e", fx)
    fx0 = @sprintf("%8.3e", fx0)
    normg = @sprintf("%7.1e", normg)
    normg0 = @sprintf("%7.1e", normg0)
    return [model.meta.name string(algo) n fx fx0 normg normg0 neval_obj(model) neval_grad(model) neval_hprod(model) k optimal] # for CUTEst problems
end
