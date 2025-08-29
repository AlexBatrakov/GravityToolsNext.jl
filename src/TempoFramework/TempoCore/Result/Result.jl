# src/TempoFramework/TempoCore/Result/Result.jl
# -----------------------------------------------------------------------------
# Aggregator for all result-related components.
# -----------------------------------------------------------------------------
include("ResidualStats.jl")
include("WhiteNoiseDiagnostics.jl")
include("WhiteNoise.jl")
include("InternalIteration.jl")
include("Convergence.jl")
include("GeneralResult.jl")