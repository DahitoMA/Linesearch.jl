using NLPModels
using LinearOperators
using OptimizationProblems
import CSV
import DataFrames

# Convex problems
# Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, brownden, cliff, clplatea, clplateb, clplatec, dixon3dq, dqdrtic, dqrtic, engval1, nasty, nondquar,
# palmer1c, palmer1d, palmer2c, palmer3c, palmer4c, palmer5c, palmer5d, palmer6c, palmer7c, palmer8c, power, quartc, tridia, vardim]

# Convex problems var ≥ 10 (17 pbs)
# Problems = [arglina, arglinb, arglinc, arwhead, bdqrtic, clplatea, clplateb, clplatec,
# dixon3dq, dqdrtic, dqrtic, engval1, nondquar, power, quartc, tridia, vardim]

# Nonconvex problems
# Problems = [AMPGO02, AMPGO03, AMPGO04, AMPGO05, AMPGO06, AMPGO08, AMPGO09,
# AMPGO10, AMPGO11, AMPGO12, AMPGO13, AMPGO14, AMPGO15, AMPGO18,
# AMPGO20, AMPGO21, AMPGO22,
# beale, brownbs, broydn7d, brybnd,
# chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
# dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
# Dus2_1, Dus2_3, Dus2_9, Duscube,
# edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
# indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
# scosine, Shpak1, Shpak2, Shpak3, Shpak4, Shpak5, Shpak6, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods];

# Nonconvex problems of dimension > 1
# Problems = [beale, brownbs, broydn7d, brybnd,
# chainwoo, chnrosnb_mod, cosine, cragglvy, curly, dixmaane, dixmaanf, dixmaang,
# dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl, dixmaanm, dixmaann, dixmaano, dixmaanp,
# edensch, eg2, errinros_mod, extrosnb, fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose, genrose_nash,
# indef_mod, liarwhd, meyer3, morebv, ncb20, ncb20b, noncvxu2, noncvxun, nondia, penalty2, powellsg, schmvett,
# scosine, sinquad, sparsine, sparsqur, srosenbr, tointgss, tquartic, woods]

# Nonconvex problems var ≥ 10 (50 pbs)
Problems = [broydn7d, brybnd, chainwoo, chnrosnb_mod, cosine, cragglvy, curly,
dixmaane, dixmaanf, dixmaang, dixmaanh, dixmaani, dixmaanj, dixmaank, dixmaanl,
dixmaanm, dixmaann, dixmaano, dixmaanp, edensch, eg2, errinros_mod, extrosnb,
fletcbv2, fletcbv3_mod, fletchcr, fminsrf2, freuroth, genhumps, genrose,
genrose_nash, indef_mod, liarwhd, morebv, ncb20, ncb20b, noncvxu2, noncvxun,
nondia, penalty2, powellsg, schmvett, scosine, sinquad, sparsine, sparsqur,
srosenbr, tointgss, tquartic, woods]


# All pbs var ≥ 10 and not in CUTEst
# Problems = [clplatea, clplateb, clplatec, chnrosnb_mod, errinros_mod,
# fletcbv3_mod, genrose_nash, scosine]

algo = CRlin
# csvname = string("noncvxlin",string(algo)[1:2],".csv")

D = ["model     algo nvar   f(x)    f(x0)   ‖g(x)‖  ‖g(x0)‖   #f  #g #Hv  #it optimal"]

MiniLogging.@info(loggerlaunch, D[1])
n = 100
for problem in Problems
    model = MathProgNLPModel(problem(n), name=string(problem))
    L = Linesearch(model, algo)
    S = @sprintf("%5s %5s %4d %5s %5s %5s %5s %4d %4d %4d %4d %5s", L[1], L[2], L[3], L[4], L[5], L[6], L[7], L[8], L[9], L[10], L[11], L[12])
    MiniLogging.@info(loggerlaunch, S)
    D = vcat(D, S)
    reset!(model)
end
# DF = DataFrames.DataFrame(D)
# CSV.write(csvname, DF)
# writedlm("ResultsCGlin.txt", D)
writedlm("nonconvpbsCRlin.txt", D)
