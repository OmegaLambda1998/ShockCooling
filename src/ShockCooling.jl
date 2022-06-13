module ShockCooling

# External packages
using Supernovae
using Unitful, UnitfulAstro
using LoggingExtras
using TOML
using Statistics
using OLUtils

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
    setup_global!(toml, verbose, Dict("base_path" => ("toml_path", ""), "output_path" => ("base_path", "Output"), "data_path" => ("base_path", "Data")))
    config = toml["global"]
    # Get data
    # Data can either be stored in the "data_path" toml file
    # Or can be specified via [ data.global ] [ data.observations ] etc...
    # TODO test
    if "data" in keys(toml)
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
    else
        supernova = nothing
    end
    # Get model
    if "model" in keys(toml)
        # You must specify the model name as it appears in src/models (i.e P15)
        @info "Loading in $(length(toml["model"])) models"
        models = [setup_model(model_dict) for model_dict in toml["model"]]
        models = Dict(model.class => model for model in models)
        @debug "Loaded in $(length(keys(models))) models"
    else
        models = nothing
    end
    if "fitting" in keys(toml)
        # Fit models to data
        priors = []
        chains = Dict()
        accept_ratios = []
        logdensities = []
        blobs = []
        if "chain" in keys(toml["fitting"])
            for model_name in keys(toml["fitting"]["chain"])
                chain_path = abspath(joinpath(config["base_path"], toml["fitting"]["chain"][model_name]))
                @info "Loading chain from $chain_path"
                chain = load_chain(chain_path)
                chains[model_name] = chain
            end
        else
            for (class, model) in models
                @info "Fitting $(model.name)"
                prior, chain, accept_ratio, logdensity, blob = run_mcmc(toml["fitting"], model, supernova)
                save_chain(joinpath(config["output_path"], "chain_$(model.name)"), chain)
                push!(priors, prior)
                chains[model.class] = chain
                push!(accept_ratios, accept_ratio)
                push!(logdensities, logdensity)
                push!(blobs, blob)
            end
        end
        bestfits = Dict()
        params = Dict()
        for (model_name, chain) in chains
            @info "Analysing $model_name"
            model = models[model_name]
            llhood = likelihood_function(model, supernova)
            bestfit = get_bestfit(chain)
            bestfits[model_name] = bestfit
            param = Dict()
            bestfit_params = [] 
            for (i, k) in enumerate(sort!(collect(keys(model.parameter_names))))
                @info "Best fitting $(model.parameter_names[k]) = $(round(bestfit[i][2], digits=3))+$(round(bestfit[i][3], digits=3))/-$(round(bestfit[i][1], digits=3)) $(model.constraints[k][2])"
                param[k] = bestfit[i][2] * model.constraints[k][2]
                push!(bestfit_params, bestfit[i][2])
            end
            likelihood = llhood(bestfit_params)
            @info "Bestfit model has likelihood $likelihood"
            params[model_name] = param
        end
    end

    if "plot" in keys(toml)
        # Plotting
        plot_config = get(toml, "plot", Dict())

        contour_plot_config = get(plot_config, "contour", nothing)
        if !isnothing(contour_plot_config)
            @info "Plotting contour"
            for (class, model) in models
                contour_plot_config["path"] = joinpath(config["output_path"], "Contour $(model.name).svg")
                plot_contour(model, chains[class], contour_plot_config) 
            end
        end

        comparison_plot_config = get(plot_config, "comparison", nothing)
        if !isnothing(comparison_plot_config)
            @info "Plotting comparison"
            for (class, model) in models
                comparison_plot_config["path"] = joinpath(config["output_path"], "Comparison $(model.name).svg")
                plot_comparison(model, supernova, params[class], chains[class], comparison_plot_config) 
            end
        end
    end
end

function run_shockcooling(toml_path::AbstractString, verbose::Bool)
    toml = TOML.parsefile(toml_path)
    if !("global" in keys(toml))
        toml["global"] = Dict()
    end
    toml["global"]["toml_path"] = dirname(abspath(toml_path))
    return run_shockcooling(toml, verbose)
end

end # module
