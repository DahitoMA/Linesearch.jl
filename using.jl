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
loggerArm = get_logger("Armijo")
loggerLin = get_logger("Linesearch")

loggerCRlin.level = MiniLogging.ERROR
loggerCGlin.level = MiniLogging.ERROR
loggerArm.level = MiniLogging.ERROR
loggerLin.level = MiniLogging.ERROR

include("CGlin.jl")
include("CRlin.jl")
include("Armijo.jl")
include("Linesearch.jl")
