using Krylov

ACG = readdlm("NewResults/AllPbsCG.txt")
ACR = readdlm("NewResults/AllPbsCR.txt")
Problems = ACG[2:end,1]

function perfprofile()
    stats = Dict{String, Any}("CR" => ["#f" "#g" "#Hv"], "CG" => ["#f" "#g" "#Hv"])
    k = 1
    for p in Problems
        k += 1
        statCR = ACR[k,12] ? ACR[k,8:10] : -ACR[k,8:10]
        stats["CR"] = vcat(values(stats["CR"]), statCR')
        statCG = ACG[k,12] ? ACG[k,8:10] : -ACG[k,8:10]
        stats["CG"] = vcat(values(stats["CG"]), statCG')
    end
    return stats
end

stats = perfprofile()
tf = font(17) # titlefont
f = font(17)
pb_type = "allpbs"
algo_used = "CR and CG"

p = performance_profile(hcat([p[2:end, 1] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="#f (all problems)",
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) #color=:black) # Profile for #f
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_f_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 2] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="#g (all problems)",
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #g
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_g_", pb_type, ".pdf"))


p = performance_profile(hcat([p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="#Hv (all problems)",
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #Hv
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
Plots.savefig(p, string("profil_Hv_", pb_type, ".pdf"))

p = performance_profile(hcat([p[2:end, 1]+p[2:end, 2]+p[2:end, 3] for p in values(stats)]...),
                            collect(String, [string(s) for s in keys(stats)]),
                            title="#f + #g + #Hv (all problems)",
                            titlefont = tf, legendfont = f, guidefont = f,
                            legend=:bottomright, palette=:blues, w=2, style=:auto) # Profile for #f + #g + #Hv
# Plots.xlabel!("Within this factor of the best (log₂ scale)")
savefig(p, string("profil_f+g+Hv_", pb_type, ".pdf"))
