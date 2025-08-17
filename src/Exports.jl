# Exports.jl
# Consolidated exports for TempoFramework

# -----------------------------
# TempoCore / AbstractTempo.jl
# -----------------------------
export AbstractTempoVersion, Tempo, Tempo2
export tempo_data_dir, tempo_cmd_path, tempo_cmd, tempo_env
export validate, get_tempo_command

# --------------------------------
# TempoCore / TempoDataFiles.jl
# --------------------------------
export TimTOAEntry, TempoResidualEntry, CombinedTOAEntry
export read_tim_file, read_residual_file, combine_tim_and_residuals

# ----------------------------
# TempoCore / TempoOutput.jl
# ----------------------------
export BasicTempoOutput, FitParameter, TempoOutputError, InternalIterationOutput
export iserror
export parse_tempo_output, parse_internal_iteration_tempo_output
export parse_basic_tempo_output, parse_fit_parameters, parse_tempo_output_error

# -------------------------------
# TempoCore / TempoParameters.jl
# -------------------------------
export GeneralTempoParameter, TP
export parse_tempo_parameter_field, extract_tempo_parameter_from_line
export get_par_file_representation
export update_or_add_tempo_parameter!, update_many_tempo_parameters!, update_params_by_dict!

# -----------------------------
# TempoCore / TempoParFile.jl
# -----------------------------
export TempoParFile
export read_par_file!, write_par_file!
export generate_par_file_name, copy_par_file
export update_one_parameter_in_par_file!, update_par_file!

# -----------------------------
# TempoCore / TempoSettings.jl
# -----------------------------
export AbstractTempoSettings
export TempoRunFiles, default_par_file_output
export TempoExecutionOptions, TempoRunModifiers, TempoRunBehavior
export WhiteNoiseAnalysisOptions
export ProcessOptions, LoggingOptions
export BasicTempoSettings
export ensure_work_dir!, validate

# -------------------------
# TempoCore / TempoRun.jl
# -------------------------
export cleanup_old_tempo_files
export run_tempo_raw, run_tempo_parsed

# ----------------------------
# TempoCore / TempoResult.jl
# ----------------------------
# basic residual stats
export BasicResidualStats, build_basic_residual_statistics
# normalized residuals
export NormalizedResidualStats, build_normalized_residual_statistics
# aggregated stats
export ResidualStatistics, build_residual_statistics
export ResidualStatisticsEntry, build_residual_statistics_entry
export ResidualStatisticsGroup, build_residual_statistics_group
export compute_global_wmean
# white-noise fit
export WhiteNoiseBackendFitResult, WhiteNoiseFitResult, build_white_noise_fit
# per-iteration result and builder
export InternalIterationResult, build_internal_iteration_result
# convergence
export ConvergenceSeries, build_convergence_series
export ConvergenceInfo, default_convergence_thresholds, is_converged, is_converged_by
export find_worst_parameter
# top-level result
export GeneralTempoResult, ParamEstimate
export build_core_metrics
export in_fit_equals_in_tim

# ---------------------------
# TempoCore / TempoTasks.jl
# ---------------------------
export AbstractTempoTask, SingleTempoTask, MultiPointTask
export run_task

# -------------------------------
# TempoCore / TempoWhiteNoise.jl
# -------------------------------
export weighted_mean, weighted_std, weighted_skewness, weighted_kurtosis, weighted_rms
export chisq_stats
export transform_uncertainties!, transform_uncertainties
export AD_objective_function!, AD_objective_function
export AD_objective_function_with_offset!, AD_objective_function_fit_offset!
export find_EFAC_for_fixed_EQUAD_and_std!
export estimate_WhiteNoise_AD_with_offset

# ---------------------------------------------
# SingleTasks / BasicTempoRun.jl (TempoFramework)
# ---------------------------------------------
export BasicTempoRun