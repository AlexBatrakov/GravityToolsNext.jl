# src/Exports.jl
# Consolidated exports for TempoFramework

# ------------------------------------------------------------------------------------------------------------------------------
# AdaptiveGridFramework
# ------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# AdaptiveGridFramework / GridAxes 
# ------------------------------------------------------------------------------------------------------------
export AbstractGridRule, LinRule, LogRule, ExplicitRule, GridAxis
export LinAxis, LogAxis, ExplicitAxis, linspace, axisvalues, refine

# ------------------------------------------------------------------------------------------------------------
# AdaptiveGridFramework / RefinementSettings
# ------------------------------------------------------------------------------------------------------------
export AbstractRefinementUnit, LocalMinimaUnit, FullUnit, DiffUnit, RelDiffUnit, ContourUnit, DiffContourUnit
export RefinementSettings

# ------------------------------------------------------------------------------------------------------------
# AdaptiveGridFramework / AdaptiveRefinement2DGrid
# ------------------------------------------------------------------------------------------------------------
export AdaptiveRefinement2DGrid, cell_selector, calculate_2DGrid, precalculate_2DGrid!, refine_2DGrid


# ------------------------------------------------------------------------------------------------------------------------------
# Utils
# ------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# Utils / PrettyPrinting 
# ------------------------------------------------------------------------------------------------------------
export print_aligned_table

# ------------------------------------------------------------------------------------------------------------
# Utils / Levels
# ------------------------------------------------------------------------------------------------------------
export lvl_1sigma, lvl_2sigma, lvl_3sigma, lvl_4sigma, lvl_5sigma, lvl_6sigma, lvl_7sigma
export lvl_68CL, lvl_90CL, lvl_95CL, lvl_99CL, lvl_997CL

# ------------------------------------------------------------------------------------------------------------------------------
# TempoFramework
# ------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------
# TempoFramework / TempoCore 
# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / AbstractTempo.jl
# ------------------------------------------------------------------------------------------
export AbstractTempoVersion,
       Tempo,
       Tempo2,
       tempo_data_dir,
       tempo_cmd_path,
       tempo_cmd,
       tempo_env,
       validate
       get_tempo_command

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoDataFiles.jl
# ------------------------------------------------------------------------------------------
export TimTOAEntry,
       TempoResidualEntry,
       CombinedTOAEntry,
       read_tim_file,
       read_tim_file_safe,
       read_residual_file,
       read_residual_file_safe,
       combine_tim_and_residuals,
       combine_tim_and_residuals_safe

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoOutput.jl
# ------------------------------------------------------------------------------------------
export BasicTempoOutput,
       FitParameter,
       TempoOutputError,
       InternalIterationOutput,
       iserror,
       parse_tempo_output
       parse_internal_iteration_tempo_output,
       parse_basic_tempo_output,
       parse_fit_parameters,
       parse_tempo_output_error

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoParameters.jl
# ------------------------------------------------------------------------------------------
export TempoParameter, TP,
       parse_tempo_parameter_field,
       extract_tempo_parameter_from_line,
       get_par_file_representation,
       value_as_float, uncertainty_as_float,
       param_index, has_param, get_param,
       upsert_param!, upsert_params!, with_upserted_params, upsert!,
       delete_param!, delete_params!, without_params

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoParFile.jl
# ------------------------------------------------------------------------------------------
export TempoParFile
export read_par_file!, read_par_file, write_par_file!
export generate_par_file_name, copy_par_file
export has_param, get_param, set_flag!
export upsert_param!, upsert_params!
export delete_param!, delete_params!, params_as_vector
export update_one_parameter_in_par_file!, update_par_file!

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoSettings.jl
# ------------------------------------------------------------------------------------------
export AbstractTempoSettings
export RunPaths, default_par_output
export EngineOptions, InputModifiers
export WhiteNoiseOptions
export WorkspaceOptions, LoggingOptions
export TempoRunSettings
export copy_with, validate, validate_inputs_exist

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoRun.jl
# ------------------------------------------------------------------------------------------
export MaterializedJob, RunArtifacts, cleanup_old_tempo_files
export run_tempo_raw, run_tempo_parsed

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / Result 
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------
# TempoFramework / TempoCore / Result / ResidualStats.jl
# ------------------------------------------------------------------------
export BasicResidualStats,
       NormalizedResidualStats,
       ResidualStatistics,
       ResidualStatisticsEntry,
       ResidualStatisticsGroup,
       build_basic_residual_statistics,
       build_normalized_residual_statistics,
       build_residual_statistics,
       build_residual_statistics_entry,
       build_residual_statistics_group,
       in_fit_equals_in_tim

# ------------------------------------------------------------------------
# TempoFramework / TempoCore / Result / WhiteNoise.jl
# ------------------------------------------------------------------------  
export WhiteNoiseBackendFitResult,
       WhiteNoiseFitResult,
       estimate_white_noise_ad_with_offset,
       build_white_noise_fit,
       print_backend_table

# ------------------------------------------------------------------------
# TempoFramework / TempoCore / Result / WhiteNoiseDiagnostics.jl
# ------------------------------------------------------------------------  
export print_white_noise_fit_report,
       print_white_noise_fit_report!

# ------------------------------------------------------------------------
# TempoFramework / TempoCore / Result / InternalIteration.jl
# ------------------------------------------------------------------------
export InternalIterationResult,
       build_internal_iteration_result

# ------------------------------------------------------------------------
# TempoFramework / TempoCore / Result / GeneralResult.jl
# ------------------------------------------------------------------------
export GeneralTempoResult,
       build_general_tempo_result,
       result_metric,
       has_converged

# ------------------------------------------------------------------------------------------
# TempoFramework / TempoCore / TempoTasks.jl
# ------------------------------------------------------------------------------------------
export AbstractTempoTask, SingleTempoTask, MultiPointTask
export run_task, task_workdir, task_with_overrides, task_with_param, task_stage_inputs!

# ------------------------------------------------------------------------------------------------------------
# TempoFramework / Prior
# ------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------
# TempoFramework / Prior / PriorSpecs.jl
# ------------------------------------------------------------------------------------------
export AbstractPriorSpec, AnalyticPrior, GridPrior, SampledPrior
export validate_prior
export prior_invcdf, prior_pdf, prior_logpdf, prior_support
export prior_median, prior_quantile, prior_quantiles
export read_prior_samples

# ------------------------------------------------------------------------------------------
# TempoFramework / Prior / NodeRules.jl
# ------------------------------------------------------------------------------------------
export AbstractNodeRule, ClenshawCurtisNodes, QuantileNodes, ExplicitThetaNodes, eval_nodes

# ------------------------------------------------------------------------------------------
# TempoFramework / Prior / PriorSettings.jl
# ------------------------------------------------------------------------------------------
export PriorExecutionOptions, PriorMarginalizationSettings, tempo_flag, build_node_dirname, write_node_metadata


# ------------------------------------------------------------------------------------------------------------
# TempoFramework / SingleTasks
# ------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# TempoFramework / SingleTasks / BasicTempoTask.jl
# -------------------------------------------------------------------------------------------
export BasicTempoTask

# -------------------------------------------------------------------------------------------
# TempoFramework / SingleTasks / PriorMarginalizedTempoTask.jl 
# -------------------------------------------------------------------------------------------
export PriorMarginalizedTempoTask


# ------------------------------------------------------------------------------------------------------------
# TempoFramework / MultiPointTasks
# ------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------
# TempoFramework / MultiPointTasks / Adaptive2DGridTask.jl
# -------------------------------------------------------------------------------------------
export Adaptive2DGridTask, GridWorkspaceOptions