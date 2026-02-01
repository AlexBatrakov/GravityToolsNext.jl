using Distributed
using GravityToolsNext
using JLD2

addprocs(8)

@everywhere using GravityToolsNext


basic_settings = TempoRunSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF",
    par_input = "DDSTG_DEF_SCINT_RN84_base.par",
    par_output = "DDSTG_DEF_SCINT_RN84_out.par",
    tim_input = "J1141-6545_pn_new.tim",
    tempo_version = Tempo2(),
    flags = "",
    nits = 5,
    gain = 1,
    override_params = [TP("STG_BETA0", 0.0), TP("STG_ALPHA0", 0.0), TP("EOS", "MPA1")],
    # time_start = nothing,
    # time_finish = nothing,
    write_output = true,
    write_residuals = true,
    save_internal_iterations = true,
    save_residuals = false,
    white_noise_enabled = true,
    white_noise_scope = :all,                     # :final, :all
    work_mode = :jobdir,                          # :inplace, :jobdir
    job_name = "MPA1_GR",                            # nothing, "TEST"
    overwrite = :clean,                           # :error, :reuse, :unique, :clean
    layout = :split,                               # :flat, :split
    temp_dir = nothing,                            # nothing, "nodes/node_1"
    link_tim = false,                              # false
    snapshot_par = true,                           # true
    cleanup_before_run = true,                     # true
    keep_tmp_on_success = false,                   # false
    keep_tmp_on_error = true,                      # true
    # timeout_s = nothing,
    write_manifest = false,                         # false
    # manifest_style = :json,                      # :json, :toml
    # verbosity = :info,                           # :silent, :warn, :info, :debug
    # with_timestamps = true,
)


basic_task = BasicTempoTask(basic_settings)

# result = run_task(basic_task)

prior_settings = PriorMarginalizationSettings(
        parameter = :DDOT,
        pin_mode = :fixed,
        prior = SampledPrior("J1141-6545_DDOT_prior.dat"),
        nodes = ClenshawCurtisNodes(4),
        likelihood_source = :chi2_fit_basic,
        ref_strategy = :prior_median,
        representative = :prior_median,
        save_node_results = true,
        exec_options = PriorExecutionOptions(
        mode               = :chained,                  # :independent | :chained
        chain_direction    = :backward,                  # :forward | :backward
        chain_snapshot_par = true,
        # scheduler          = :serial,                   # :serial | :threads | :distributed
        # max_workers        = 0,
        # workdir_layout     = :per_node,                 # :per_node
        # node_dir_prefix    = "nodes/node_",
        # keep_node_dirs     = true,
        on_error             = :stop,
        # dir_name_mode      = :index_only,               # :index_only | :with_value
        # index_pad          = 3,
        # value_sig          = 6
        )
)

prior_task = PriorMarginalizedTempoTask(basic_task, prior_settings)

# prior_result = run_task(prior_task)
# jldsave("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/MPA1_GR/prior_result.jld2"; prior_result = prior_result)
prior_result = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/MPA1_GR/prior_result.jld2", "prior_result")

prior_settings = copy_with(prior_settings; 
                            seed_spec = SeedPaths(prior_result),
                            exec_options = copy_with(prior_settings.exec_options; mode = :independent))

prior_task = PriorMarginalizedTempoTask(basic_task, prior_settings)

grid_task = Adaptive2DGridTask(
    base_task = prior_task,
    x = GridAxis(:STG_BETA0,  min=-6.0, max=6.0,    N=4, rule = LinRule()),     # name, min, max, initial nodes
    y = GridAxis(:STG_ALPHA0, min=-1e-4, max=-1e-1, N=4, rule = LogRule()),     # name, min, max, initial nodes
    ref_settings = RefinementSettings(
        params_to_save = (:chi2_marginalized, :wrms_fit, :wrms_tn_fit, :ad_white_fit, :pre_post_final, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT),
        desired_refinement_level = 0,
        parallel = true,
        DiffContourUnit(:chi2_marginalized, from_min=true, diffs=[1.0], contours=[lvl_3sigma], nan_default=Inf),
        LocalMinimaUnit(:chi2_marginalized, from_min=true, max=20.0, max_diff = 0.1, nan_default=Inf)
    ),
    opts = GridWorkspaceOptions(
        grid_root = "MPA1_GRID"
    ),
)

# jldsave("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/FAKE_TEST/grid_result.jld2"; grid_result = grid_result)
# grid_result = load("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/FAKE_TEST/grid_result.jld2", "grid_result")

base_path = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_DEF/MPA1_GRID"
# grid_result = load(joinpath(base_path, "grid_result_iter3.jld2"), "grid_result")


grid_result = run_task(grid_task)
jldsave(joinpath(base_path, "grid_result_iter0.jld2"); grid_result = grid_result)
plot_chi2_contours(grid_result)

for i in 3:3
    grid_result = run_task(grid_task; grid_init = grid_result, just_refine = true)
    jldsave(joinpath(base_path, "grid_result_iter$(i).jld2"); grid_result = grid_result)
    plot_chi2_contours(grid_result)
end

using PyPlot
using LaTeXStrings

function plot_chi2_contours(grid_result::AdaptiveRefinement2DGrid, param_name::Symbol = :chi2_marginalized)

    chi2_min = grid_result.min[param_name]
    rc("mathtext",fontset="cm")
    rc("font", family="serif", size=12)
    # y_log10_values = log10.(-grid_result.y.values)
    fig, ax = subplots()
    y_delta = grid_result.y.lin_values[2] - grid_result.y.lin_values[1]
    x_delta = grid_result.x.values[2] - grid_result.x.values[1]
    extent=(grid_result.x.values[1] - 0.5*x_delta, grid_result.x.values[end] + 0.5*x_delta, grid_result.y.lin_values[1] - 0.5*y_delta, grid_result.y.lin_values[end] + 0.5*y_delta)
    imsh = ax.imshow(permutedims(grid_result.vars[param_name] .- grid_result.min[param_name]), extent=extent, cmap="Blues_r", norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20.0), origin="lower", aspect="auto")
    cbar = colorbar(imsh)
    cs = ax.contour(grid_result.x.values, grid_result.y.lin_values, permutedims(grid_result.vars[param_name] .- grid_result.min[param_name]), levels=[lvl_1sigma, lvl_2sigma, lvl_3sigma], linestyles=["-", "--", "-."], colors="red")
    plot([], [], label=L"\Delta\chi^{2} (1\sigma, 2\sigma, 3\sigma)", "-",  color="red")
    chi2_min_round = round(chi2_min, digits=2)
    annot_text = L"$\chi^2_\mathrm{min} = %$chi2_min_round$"
    # Добавление текста с LaTeX
    text(
        0.05, 0.05, annot_text, transform=gca().transAxes,
        fontsize=12,
        bbox=Dict("facecolor" => "white", "alpha" => 0.8, "boxstyle" => "round")
    )
    ax.set_ylabel(L"\log_{10}|\alpha_0|", size=14)
    ax.set_xlabel(L"\beta_0", size=14)
    # title("Red noise + white noise")
    cbar.set_label(L"$\Delta\chi^{2} \equiv \chi^{2} - \chi^2_\mathrm{min}$", fontsize=14)
    legend(fontsize=12)
    tight_layout()

end