module ShockCooling

# External packages
using Supernovae
using Unitful, UnitfulAstro
using LoggingExtras
using TOML
using Statistics

# Internal Packages
include("Models.jl")
using .Models

include("Fitting.jl")
using .Fitting

include("Plotting.jl")
using .Plotting

# Exports
export run_shockcooling
export plot_luminosity, plot_luminosity!
export plot_radius, plot_radius!
export plot_temperature, plot_temperature!
export plot_prior, plot_prior!
export plot_walkers, plot_walkers!
export plot_contour, plot_contour!
export plot_comparison, plot_comparison!

function setup_global_config!(toml::Dict)
    config = get(toml, "global", Dict())
    # Base path is where everything relative will be relative to
    # Defaults to the directory containing the toml path
    # Can be relative (to the toml path) or absolute
    base_path = get(config, "base_path", nothing)
    if isnothing(base_path)
        base_path = dirname(toml["toml_path"])
    elseif !isabspath(base_path)
        base_path = joinpath(dirname(toml["toml_path"]), base_path)
    end
    base_path = abspath(base_path)
    config["base_path"] = base_path
    # Output path is where all output (figures) will be placed
    # Defaults to base_path / Output
    # Can be relative (to base_path) or absolute
    output_path = get(config, "output_path", nothing)
    if isnothing(output_path)
        output_path = joinpath(base_path, "Output")
    elseif !isabspath(output_path)
        output_path = joinpath(base_path, output_path)
    end
    config["output_path"] = abspath(output_path)
    # Data path is where all supernovae data (both photometric and spectroscopic) will be stored
    # Default to base_path / Data
    # Can be relatvie (to base_path) or absolute
    data_path = get(config, "data_path", nothing)
    if isnothing(data_path)
        data_path = joinpath(base_path, "Data")
    elseif !isabspath(data_path)
        data_path = joinpath(base_path, data_path)
    end
    config["data_path"] = abspath(data_path)
    # Logging sets whether or not to setup and use Supernovae's logging
    logging = get(config, "logging", false)
    config["logging"] = logging
    # Log file is the name of the log file. This will only work if logging is true
    # Can only be relative to output_path
    # Defaults to log.txt
    log_file = get(config, "log_file", nothing)
    if logging
        if isnothing(log_file)
            log_file = "log.txt"
        end
        log_file = abspath(joinpath(output_path, log_file))
    end
    if !logging & !isnothing(log_file)
        @warn "Logging set to false, so log file $log_file will not be written. Please add `logger=true` to your [ global ] config"
    end
    config["log_file"] = log_file
    toml["global"] = config
end

function setup_logger(log_file::AbstractString, verbose::Bool)
    if verbose
        level = Logging.Debug
    else
        level = Logging.Info
    end
    function fmt(io, args)
        if args.level == Logging.Error
            color = :red
            bold = true
        elseif args.level == Logging.Warn
            color = :yellow
            bold = true
        elseif args.level == Logging.Info
            color = :cyan
            bold = false
        else
            color = :white
            bold = false
        end
        printstyled(io, args._module, " | ", "[", args.level, "] ", args.message, "\n"; color = color, bold = bold)
    end
    logger = TeeLogger(
        MinLevelLogger(FormatLogger(fmt, open(log_file, "w")), level),
        MinLevelLogger(FormatLogger(fmt, stdout), level)
    )
    global_logger(logger)
    @info "Logging to $log_file"
end

function setup_model(model_dict::Dict)
    # TODO Add the ability to change parameter names
    model = get_model(model_dict["name"])
    constraints = get_constraints(model_dict["constraints"])
    return model(constraints)
end

function update_supernova!(supernova::Supernova, config::Dict)
    @debug "Starting with $(length(supernova.lightcurve.observations)) observations"
    # Unpack
    min_time = get(config, "min_time", -Inf)
    min_time_unit = uparse(get(config, "min_time_unit", "d"), unit_context = [Unitful, UnitfulAstro])
    min_time_range = get(config, "min_time_range", [-Inf, Inf]) .* min_time_unit
    max_time = get(config, "max_time", Inf)
    max_time_unit = uparse(get(config, "max_time_unit", "d"), unit_context = [Unitful, UnitfulAstro])
    max_time_range = get(config, "max_time_range", [-Inf, Inf]) .* max_time_unit

    # Get min time
    @debug "min_time option = $min_time"
    if min_time in ["min", "max"]
        observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > min_time_range[1]) & (obs.time < min_time_range[2]))]
        best_obs = observations[1]
        for obs in observations
            if ((min_time == "min") & (obs.flux < best_obs.flux)) | ((min_time == "max") & (obs.flux > best_obs.flux))
                best_obs = obs
            end
        end
        min_time = best_obs.time
    else
        min_time = min_time * min_time_unit
    end
    @debug "Min time set to $min_time"

    # Get max time
    @debug "max_time option = $max_time"
    if max_time in ["min", "max"]
        observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > max_time_range[1]) & (obs.time < max_time_range[2]))]
        best_obs = observations[1]
        for obs in observations
            if ((max_time == "min") & (obs.flux < best_obs.flux)) | ((max_time == "max") & (obs.flux > best_obs.flux))
                best_obs = obs
            end
        end
        max_time = best_obs.time
    else
        max_time = max_time * max_time_unit
    end
    @debug "Max time set to $max_time"

    # Cut out observations
    observations = [obs for obs in supernova.lightcurve.observations if ((obs.time > min_time) & (obs.time < max_time))]
    supernova.lightcurve.observations = observations
    @debug "Ending with $(length(supernova.lightcurve.observations)) observations"
end

function get_bestfit(chain)
    m, r, t, v = [quantile([c[i] for c in chain], [0.16, 0.5, 0.84]) for i in 1:4]
    for val in [m, r, t, v]
        val[1] = val[2] - val[1]
        val[3] = val[3] - val[2]
    end
    return m, r, t, v
end

function run_shockcooling(toml::Dict, verbose::Bool)
    setup_global_config!(toml)
    config = toml["global"]
    # Ensure all path's exist
    # TODO test
    if !isdir(config["base_path"])
        mkpath(config["base_path"])
    end
    if !isdir(config["output_path"])
        mkpath(config["output_path"])
    end
    # Optionally set up logging
    if config["logging"]
        setup_logger(config["log_file"], verbose)
    end
    # Get data
    # Data can either be stored in the "data_path" toml file
    # Or can be specified via [ data.global ] [ data.observations ] etc...
    # TODO test
    @info "Loading in data"
    data_path = get(toml["data"], "data_path", nothing)
    if isnothing(data_path)
        @debug "Using data in original toml"
        # Assuming all relevant info is in the [ data ] block
        supernova = process_supernova(toml["data"], verbose)
    else
        if !isabspath(data_path)
            data_path = joinpath(config["base_path"], data_path)
        end
        @debug "Using data in $data_path"
        supernova = process_supernova(data_path, verbose)
    end
    # Cut down data
    update_supernova!(supernova, toml["data"])
    # Get model
    # You must specify the model name as it appears in src/models (i.e P15)
    @info "Loading in $(length(toml["model"])) models"
    models = [setup_model(model_dict) for model_dict in toml["model"]]
    @debug "Loaded in $(length(models)) models"
    # Fit models to data
    priors = []
    chains = []
    accept_ratios = []
    logdensities = []
    blobs = []
    bestfits = []
    params = []
    for model in models
        @info "Fitting $(model.name)"
        prior, chain, accept_ratio, logdensity, blob = run_mcmc(toml["fitting"], model, supernova)
        save_chain(joinpath(config["output_path"], "chain_P15"), chain)
        push!(priors, prior)
        push!(chains, chain)
        push!(accept_ratios, accept_ratio)
        push!(logdensities, logdensity)
        push!(blobs, blob)
        bestfit = get_bestfit(chain)
        push!(bestfits, bestfit)
        param = Dict()
        for (i, k) in enumerate(sort!(collect(keys(model.parameter_names))))
            @info "Best fitting $(model.parameter_names[k]) = $(bestfit[i][2])+$(bestfit[i][3])/-$(bestfit[i][1]) $(model.constraints[k][2])"
            param[k] = bestfit[i][2] * model.constraints[k][2]
        end
        push!(params, param)
    end

    # Plotting
    plot_config = get(toml, "plot", Dict())

    contour_plot_config = get(plot_config, "contour", nothing)
    if !isnothing(contour_plot_config)
        @info "Plotting contour"
        for (i, model) in enumerate(models)
            contour_plot_config["path"] = joinpath(config["output_path"], "contour_$(model.name).svg")
            plot_contour(model, chains[i], contour_plot_config) 
        end
    end

    comparison_plot_config = get(plot_config, "comparison", nothing)
    if !isnothing(comparison_plot_config)
        @info "Plotting comparison"
        for (i, model) in enumerate(models)
            comparison_plot_config["path"] = joinpath(config["output_path"], "comparison_$(model.name).svg")
            plot_comparison(model, supernova, params[i], comparison_plot_config) 
        end
    end
end

function run_shockcooling(toml_path::AbstractString, verbose::Bool)
    toml = TOML.parsefile(toml_path)
    toml["toml_path"] = abspath(toml_path)
    return run_shockcooling(toml, verbose)
end

end # module
