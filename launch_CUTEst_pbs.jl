using NLPModels
using CUTEst
include("Armijo.jl")
include("LinesearchCUTEst.jl")
include("CGlin.jl")
include("CRlin.jl")

num = parse(Int64,ARGS[1])
Problems = open(readlines, "Problems.txt")
# Problems = deleteat!(Problems,56) # pb "INDEF" not treated
Algos = [CGlin, CRlin]

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g #Hv  #it optimal"]

problem = Problems[num]
println(problem)
model = CUTEstModel(problem)
for algo in Algos
    L = LinesearchCUTEst(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %5s", L[1], L[2], L[3], L[4], L[5], L[6], L[7], L[8], L[9], L[10], L[11], L[12])
    D = vcat(D, S)
    reset!(model)
end
finalize(model)
writedlm(string(problem, "_CUTEst.txt"), D)
println("resolution ok !!!")
