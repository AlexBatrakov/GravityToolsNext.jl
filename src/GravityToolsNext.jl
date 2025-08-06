module GravityToolsNext

# --- Базовые пакеты ---
using Printf
using DelimitedFiles
using StructArrays
using Statistics
using StatsBase
using Distributions
using HypothesisTests

# using LinearAlgebra
# using Random

# using Distributed
# using PyPlot
# using ProgressMeter
# using QuadGK


# Подмодули ядра
include("TempoFramework/TempoCore/AbstractTempo.jl")
include("TempoFramework/TempoCore/TempoParameters.jl")
include("TempoFramework/TempoCore/TempoParFile.jl")
include("TempoFramework/TempoCore/TempoDataFiles.jl")
include("TempoFramework/TempoCore/TempoOutput.jl")
include("TempoFramework/TempoCore/TempoResult.jl")
include("TempoFramework/TempoCore/TempoTasks.jl")
include("TempoFramework/TempoCore/TempoSettings.jl")
include("TempoFramework/TempoCore/TempoRun.jl")

# Настройки
# include("TempoFramework/Settings/BasicTempoSettings.jl")

# Запуск
include("TempoFramework/SingleTasks/BasicTempoRun.jl")

include("Exports.jl")

end # module GravityToolsNext
