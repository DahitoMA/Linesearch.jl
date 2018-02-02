using NLPModels
using LinearOperators
using Plots # graph
using BenchmarkTools
using BenchmarkProfiles
using MiniLogging

# basic configuration. The root logger level is then INFO
basic_config(MiniLogging.INFO; date_format="%Y-%m-%d %H:%M:%S")
loggerCRlin = get_logger("CRlin")
loggerCGlin = get_logger("CGlin")
loggerWolfe = get_logger("StrongWolfe")
loggerLin = get_logger("Linesearch")
loggerlaunch = get_logger("launch_linesearch.jl")

loggerCRlin.level = MiniLogging.ERROR
loggerCGlin.level = MiniLogging.ERROR
loggerWolfe.level = MiniLogging.ERROR
loggerLin.level = MiniLogging.ERROR
# loggerlaunch.level = MiniLogging.ERROR

include("CGlin.jl")
include("CRlin.jl")
include("StrongWolfe.jl")
include("Linesearch.jl")
include("LCG.jl")
include("LCR.jl")
