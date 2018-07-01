using OptimizationProblems

## All pbs dimension > 1 and not in CUTEst
Problems = [brownden, cliff, clplatea, clplateb, clplatec, nasty,
palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d,
palmer6c, palmer7c, palmer8c,
beale, brownbs,
chnrosnb_mod,
errinros_mod, fletcbv3_mod, genrose_nash,
meyer3, scosine]

algo = CGlin
# algo = CRlin

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g  #Hv  #it #s_it #vs_it #rej_it optimal"]

for problem in Problems
    model = MathProgNLPModel(problem(), name=string(problem))
    L = Linesearch(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %5s", L[1], L[2], L[3], L[4], L[5], L[6], L[7], L[8], L[9], L[10], L[11], L[12])
    D = vcat(D, S)
    reset!(model)
end
writedlm("PbsCG.txt", D)
