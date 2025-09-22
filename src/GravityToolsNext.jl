module GravityToolsNext

# --- Базовые пакеты ---
using Printf
using Unicode
using DelimitedFiles
using StructArrays
using Statistics
using StatsBase
using Distributions
using KernelDensity
using HypothesisTests
using Roots
using Optim
using Dierckx
using QuadGK
using Dates
using Random
using JLD2
using SHA

using Distributed
# using PyPlot
using ProgressMeter
# using QuadGK

# Utils
include("Utils/PrettyPrinting.jl")

# AdaptiveGridFramework
include("AdaptiveGridFramework/GridAxes.jl")
include("AdaptiveGridFramework/RefinementSettings.jl")
include("AdaptiveGridFramework/AdaptiveRefinement2DGrid.jl")

# TempoFramework/TempoCore
include("TempoFramework/TempoCore/AbstractTempo.jl")
include("TempoFramework/TempoCore/TempoParameters.jl")
include("TempoFramework/TempoCore/TempoParFile.jl")
include("TempoFramework/TempoCore/TempoDataFiles.jl")
include("TempoFramework/TempoCore/TempoOutput.jl")
include("TempoFramework/TempoCore/Result/Result.jl")
include("TempoFramework/TempoCore/TempoTasks.jl")
include("TempoFramework/TempoCore/TempoSettings.jl")
include("TempoFramework/TempoCore/TempoRun.jl")

# TempoFramework/Prior
include("TempoFramework/Prior/PriorSpecs.jl")
include("TempoFramework/Prior/NodeRules.jl")
include("TempoFramework/Prior/PriorSettings.jl")

# TempoFramework/SingleTasks
include("TempoFramework/SingleTasks/BasicTempoTask.jl")
include("TempoFramework/SingleTasks/IterativeTempoTask.jl")
include("TempoFramework/SingleTasks/PriorMarginalizedTempoTask.jl")

# TempoFramework/MultiPointTasks
include("TempoFramework/MultiPointTasks/Adaptive2DGridTask.jl")

include("Exports.jl")

end # module GravityToolsNext
