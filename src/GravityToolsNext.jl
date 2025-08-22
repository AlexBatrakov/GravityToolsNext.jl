module GravityToolsNext

# --- Базовые пакеты ---
using Printf
using DelimitedFiles
using StructArrays
using Statistics
using StatsBase
using Distributions
using KernelDensity
using HypothesisTests
using LinearAlgebra
using Roots
using Optim
# using Random

# using Distributed
# using PyPlot
# using ProgressMeter
# using QuadGK


# TempoFramework/TempoCore
include("TempoFramework/TempoCore/AbstractTempo.jl")
include("TempoFramework/TempoCore/TempoParameters.jl")
include("TempoFramework/TempoCore/TempoParFile.jl")
include("TempoFramework/TempoCore/TempoDataFiles.jl")
include("TempoFramework/TempoCore/TempoOutput.jl")
include("TempoFramework/TempoCore/TempoResult.jl")
include("TempoFramework/TempoCore/TempoTasks.jl")
include("TempoFramework/TempoCore/TempoSettings.jl")
include("TempoFramework/TempoCore/TempoRun.jl")
include("TempoFramework/TempoCore/TempoWhiteNoise.jl")

# TempoFramework/Prior
include("TempoFramework/Prior/PriorSpecs.jl")
include("TempoFramework/Prior/NodeRules.jl")
include("TempoFramework/Prior/PriorSettings.jl")

# TempoFramework/SingleTasks
include("TempoFramework/SingleTasks/BasicTempoTask.jl")
include("TempoFramework/SingleTasks/IterativeTempoTask.jl")
include("TempoFramework/SingleTasks/PriorMarginalizedTempoTask.jl")

include("Exports.jl")

end # module GravityToolsNext
