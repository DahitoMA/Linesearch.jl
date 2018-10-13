using NLPModels
using LinearOperators
using Plots # graph
Plots.pyplot()
using BenchmarkTools
using BenchmarkProfiles
using MiniLogging

# basic configuration. The root logger level is then INFO
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")
# loggerCRlin = MiniLogging.Logger("CRlin", "CRlin.out", MiniLogging.DEBUG)
# loggerCRlin = MiniLogging.Logger("CRlin", "CRlin.out", MiniLogging.DEBUG)

loggerCRlin = get_logger("CRlin")
loggerCGlin = get_logger("CGlin")
loggerArm = get_logger("Armijo")
loggerLin = get_logger("Linesearch")
loggerlaunch = get_logger("launch_linesearch.jl")
loggerLinCUTEst = get_logger("LinesearchCUTEst")

loggerCRlin.level = MiniLogging.ERROR
loggerCGlin.level = MiniLogging.ERROR
loggerArm.level = MiniLogging.ERROR
loggerLin.level = MiniLogging.ERROR
# loggerlaunch.level = MiniLogging.ERROR
loggerLinCUTEst.level = MiniLogging.ERROR


include("CGlin.jl")
include("CRlin.jl")
include("Armijo.jl")
include("Linesearch.jl")
include("LCG.jl")
include("LCR.jl")
include("LinesearchCUTEst.jl")
