#------------------------------------------------------------------------------------------------------------------------------------

function weighted_mean(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
	weights = unc .^ -2
	sum_weights = sum(weights)
	weighted_mean_res = sum(weights .* res) / sum_weights
	weighted_mean_unc = sqrt(1 / sum_weights)
	return weighted_mean_res, weighted_mean_unc
end

function weighted_std(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
	weights = unc .^ -2
	sum_weights = sum(weights)
	weighted_mean_res = sum(weights .* res) / sum_weights
	weighted_std_res = sqrt(sum(weights .* (res .- weighted_mean_res).^2) / sum_weights)
	N = length(res)
	N_eff = (sum_weights^2) / sum(weights.^2)
	weighted_std_res_unc = weighted_std_res / sqrt(2 * (N_eff - 1))
	return weighted_std_res, weighted_std_res_unc
end

function weighted_skewness(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)
    N = length(res)

    # Взвешенное среднее
    mean_w = sum(weights .* res) / sum_weights

    # Взвешенное стандартное отклонение
    std_w = sqrt(sum(weights .* (res .- mean_w).^2) / sum_weights)

    # Взвешенный skewness
    skew_numer = sum(weights .* (res .- mean_w).^3)
    skew_denom = sum_weights * std_w^3
    skew_w = skew_numer / skew_denom

    # Эффективное число наблюдений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка skewness
    skew_err = sqrt(6.0 / N_eff)

    return skew_w, skew_err
end

function weighted_kurtosis(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)
    N = length(res)

    # Взвешенное среднее
    mean_w = sum(weights .* res) / sum_weights

    # Взвешенное стандартное отклонение
    std_w = sqrt(sum(weights .* (res .- mean_w).^2) / sum_weights)

    # Взвешенный kurtosis (не эксцесс)
    kurt_numer = sum(weights .* (res .- mean_w).^4)
    kurt_denom = sum_weights * std_w^4
    kurt_w = kurt_numer / kurt_denom

    # Эффективное число наблюдений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка kurtosis
    kurt_err = sqrt(24.0 / N_eff)

    # Возвращаем также эксцесс
    excess_kurt_w = kurt_w - 3.0

    return excess_kurt_w, kurt_err
end

function weighted_rms(res::Vector{Float64}, unc::Vector{Float64} = ones(length(res)))
    weights = unc .^ -2
    sum_weights = sum(weights)

    # RMS от нуля
    rms = sqrt(sum(weights .* res.^2) / sum_weights)

    # Эффективное число измерений
    N_eff = sum_weights^2 / sum(weights.^2)

    # Ошибка RMS
    rms_err = rms / sqrt(2 * N_eff)

    return rms, rms_err
end

function chisq_stats(res::Vector{Float64}, unc::Vector{Float64}, dof::Int = length(res))
    # Предполагается, что unc уже включает EFAC и EQUAD
    chi2 = sum((res ./ unc).^2)
    red_chi2 = chi2 / dof
    return chi2, red_chi2
end

function transform_uncertainties!(uncertainties_transformed::Vector{Float64}, uncertainties_orig::Vector{Float64}, EFAC::Float64, EQUAD::Float64)
    @. uncertainties_transformed = sqrt(EFAC^2 * uncertainties_orig^2 + EQUAD^2)
    return uncertainties_transformed
end

function transform_uncertainties(uncertainties_orig::Vector{Float64}, EFAC::Float64, EQUAD::Float64)
	uncertainties_transformed = similar(uncertainties_orig)
	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
	return uncertainties_transformed
end


function AD_objective_function!(
    residuals_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64}
	)
    @. residuals_norm = residuals / uncertainties_transformed
    ad_test = OneSampleADTest(residuals_norm, Normal(0, 1))
    return ad_test.A²
end

function AD_objective_function(residuals_norm::Vector{Float64})
    ad_test = OneSampleADTest(residuals_norm, Normal(0, 1))
    return ad_test.A²
end

function AD_objective_function_with_offset!(
    residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64},
    offset::Float64
	)
    @. residuals_shifted_norm = (residuals - offset) / uncertainties_transformed
    ad_test = OneSampleADTest(residuals_shifted_norm, Normal(0, 1))
    return ad_test.A²
end


function AD_objective_function_fit_offset!(
	residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_transformed::Vector{Float64}
	)

    function AD_objective_local(args)
    	offset = args[1]
        return AD_objective_function_with_offset!(
            residuals_shifted_norm,
            residuals,
            uncertainties_transformed,
            offset
        )
    end

    optim_result = optimize(AD_objective_local, [0.0])
    offset = optim_result.minimizer[1]
    AD_objective_val = AD_objective_local(offset)
    return AD_objective_val, offset
end

function find_EFAC_for_fixed_EQUAD_and_std!(
	residuals_norm::Vector{Float64},
	uncertainties_transformed::Vector{Float64},
	residuals::Vector{Float64},
	uncertainties_orig::Vector{Float64},
	EQUAD_fixed::Float64,
	std_res_norm_fixed::Float64
	)

    function std_res_norm_objective(EFAC::Float64)
        transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD_fixed)
        @. residuals_norm = residuals / uncertainties_transformed
        return std_res_norm_fixed / std(residuals_norm) - 1
    end

    std_residuals = std(residuals)

    if EQUAD_fixed > std_residuals / std_res_norm_fixed
        return NaN
    end

    mean_uncertainties_orig_squared = dot(uncertainties_orig, uncertainties_orig) / length(uncertainties_orig)

    EFAC_init = sqrt(((std_residuals / std_res_norm_fixed)^2 - EQUAD_fixed^2) / mean_uncertainties_orig_squared)

    EFAC = NaN

    try 
    	EFAC = find_zero(std_res_norm_objective, EFAC_init)
    catch error 
    	println(residuals)
    	println(uncertainties_orig)
    	println(EQUAD_fixed)
    	println(std_res_norm_fixed)
    	readline()
    end
    return abs(EFAC)
end

function estimate_WhiteNoise_AD_with_offset(residuals, uncertainties_orig; print_results = false, plot_results = false, backend = "", EFAC_init_in = nothing, EQUAD_init_in = nothing, offset_init_in = nothing)
	N_TOAs = length(residuals)
	SD_mean = 1.0 / sqrt(N_TOAs)
	SD_median = 1.0 / sqrt(2.0 * pi * N_TOAs)
	SD_std = 1.0 / sqrt(2 * (N_TOAs - 1))
	SD_skewness = sqrt(6.0 / N_TOAs)
	SD_kurtosis = sqrt(24.0 / N_TOAs)

	EFAC_max  = std(residuals ./ uncertainties_orig)
	EQUAD_max = std(residuals)

	uncertainties_transformed = similar(uncertainties_orig)
	residuals_norm            = similar(residuals)
	residuals_shifted_norm    = similar(residuals)
	
	function make_objective_function(
			residuals_shifted_norm::Vector{Float64},
			uncertainties_transformed::Vector{Float64},
			residuals::Vector{Float64},
			uncertainties_orig::Vector{Float64},
			)

    	function objective_fun(args)
        	EFAC, EQUAD, offset = args

        	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)

        	AD_objective = AD_objective_function_with_offset!(
            	residuals_shifted_norm,
            	residuals,
            	uncertainties_transformed,
            	offset
        	)

        	return AD_objective
    	end

    	return objective_fun
	end

	EFAC_EQUAD_offset_AD_objective = make_objective_function(residuals_shifted_norm, uncertainties_transformed, residuals, uncertainties_orig)


    if isnothing(EFAC_init_in) || isnothing(EQUAD_init_in) || isnothing(offset_init_in)
    	N_EQUADs = 11
    	N_std_res_norm = 11

    	EQUAD_init_arr = collect(LinRange(0.0, EQUAD_max, N_EQUADs))
    	std_res_norm_init_arr   = collect(LinRange(0.95, 1.05, N_std_res_norm))

    	EFAC_init_best, EQUAD_init_best, offset_init_best, AD_objective_init_best = NaN, NaN, NaN, Inf

    	for EQUAD_init in EQUAD_init_arr, std_res_norm_init in std_res_norm_init_arr 

    		EFAC_init = find_EFAC_for_fixed_EQUAD_and_std!(residuals_norm, uncertainties_transformed, residuals, uncertainties_orig, EQUAD_init, std_res_norm_init)

    		if isnan(EFAC_init)
    			continue
    		end

    		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC_init, EQUAD_init)

    		AD_objective_init, offset_init = AD_objective_function_fit_offset!(residuals_shifted_norm, residuals, uncertainties_transformed)

    		if AD_objective_init < AD_objective_init_best
    			AD_objective_init_best = AD_objective_init
    			EFAC_init_best         = EFAC_init
    			EQUAD_init_best        = EQUAD_init
    			offset_init_best       = offset_init
    		end
    	end

    	AD_result = optimize(EFAC_EQUAD_offset_AD_objective, [EFAC_init_best, EQUAD_init_best, offset_init_best])

    else
    	AD_result = optimize(EFAC_EQUAD_offset_AD_objective, [EFAC_init_in, EQUAD_init_in, offset_init_in])
    end

    EFAC_best, EQUAD_best, offset_best = Optim.minimizer(AD_result)
    EFAC_best  = abs(EFAC_best)
    EQUAD_best = abs(EQUAD_best)

    AD_objective_best = AD_result.minimum

    if print_results

    	transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC_best, EQUAD_best)
    	@. residuals_norm = residuals / uncertainties_transformed

    	println("backend: $backend")
    	println("	N_TOAs: $N_TOAs")
		println("residuals stats:")
		println("	weighted mean:     $(weighted_mean(residuals, uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(residuals, uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(residuals, uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(residuals, uncertainties_transformed))")

		println("original uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_orig))")
		println("	weighted std:      $(weighted_std(uncertainties_orig))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_orig))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_orig))")

		println("transformed uncertainties stats:")
		println("	weighted mean:     $(weighted_mean(uncertainties_transformed))")
		println("	weighted std:      $(weighted_std(uncertainties_transformed))")
		println("	weighted skewness: $(weighted_skewness(uncertainties_transformed))")
		println("	weighted kurtosis: $(weighted_kurtosis(uncertainties_transformed))")

		println("normalized residuals stats:")
		println("	mean:     $(weighted_mean(residuals_norm))")
		println("	std:      $(weighted_std(residuals_norm))")
		println("	skewness: $(weighted_skewness(residuals_norm))")
		println("	kurtosis: $(weighted_kurtosis(residuals_norm))")

		println("chi squared stats:")
		chi2, chi2r = chisq_stats(residuals, uncertainties_transformed)
		println("	chi2: $chi2")
		println("	chi2r: $chi2r")

		println("AD_objective_best:  $AD_objective_best")
		println("EFAC_best:   $EFAC_best")
		println("EQUAD_best:  $EQUAD_best")
		println("offset_best: $offset_best")

		println()

    end

    if plot_results
    	
		EFAC_arr  = collect(LinRange(0.01, ceil(EFAC_max), 128))
		EQUAD_arr = collect(LinRange(0.0, ceil(EQUAD_max), 128))
		AD_objective_arr = [
    		begin
        		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
        		AD_objective_function_fit_offset!(residuals_shifted_norm, residuals, uncertainties_transformed)[1]
    		end
    		for EFAC in EFAC_arr, EQUAD in EQUAD_arr
		]

		std_residuals_norm_arr = [
    		begin
        		transform_uncertainties!(uncertainties_transformed, uncertainties_orig, EFAC, EQUAD)
        		@. residuals_norm = residuals / uncertainties_transformed
        		std(residuals_norm)
    		end
    		for EFAC in EFAC_arr, EQUAD in EQUAD_arr
		]

		fig, ax = subplots()
		imsh = ax.imshow(log10.(AD_objective_arr), extent=(EQUAD_arr[1], EQUAD_arr[end], EFAC_arr[1], EFAC_arr[end]), cmap="Blues_r", origin="lower", aspect="auto")
		cbar = colorbar(imsh)
		cbar.set_label(L"$\log_{10}\mathrm{AD_objective}$", fontsize=14)

		# cs2 = ax.contour(EQUAD_arr, EFAC_arr, AD_objective_arr, levels = [0.5, 0.75, 1.0, 1.25], linestyles=["-", "--", "-.", ":"], colors="green")
		# plot([], [], label="AD statistics", "-",  color="green")

		cs2 = ax.contour(EQUAD_arr, EFAC_arr, AD_objective_arr, levels = AD_objective_best .+ [0.1, 0.2, 0.3, 0.4], linestyles=["-", "--", "-.", ":"], colors="green")
		plot([], [], label="AD statistics", "-",  color="green")

		cs3 = ax.contour(EQUAD_arr, EFAC_arr, std_residuals_norm_arr, levels = [1.0], linestyles=["-"], colors="violet")
		plot([], [], label="std", "-",  color="violet")

		plot(EQUAD_best, EFAC_best, "x", color="red")

		AD_objective_grid, min_ind = findmin(AD_objective_arr)
		plot(EQUAD_arr[min_ind[2]], EFAC_arr[min_ind[1]], "x", color="black")

		println("AD_objective: optimized = $AD_objective_best, grid = $AD_objective_grid, delta = $(AD_objective_best - AD_objective_grid)")
		println("std: optimized = $(std(residuals_norm)), grid = $(std_residuals_norm_arr[min_ind])")

		xlabel("EQUAD")
		ylabel("EFAC")
		title("$backend")
		legend()
		tight_layout()

	end

    return EFAC_best, EQUAD_best, offset_best, AD_objective_best
end


