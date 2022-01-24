module ShockCooling

# External packages
using Supernovae
using LoggingExtras
using TOML

# Internal Packages
include("Models.jl")
using .Models

# Exports
export run_shockcooling

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
    # Get model
    # You must specify the model name as it appears in src/models (i.e P15)
    @info "Loading in $(length(toml["model"])) models"
    models = [setup_model(model_dict) for model_dict in toml["model"]]
    @debug "Loaded in $(length(models)) models"
end

function run_shockcooling(toml_path::AbstractString, verbose::Bool)
    toml = TOML.parsefile(toml_path)
    toml["toml_path"] = abspath(toml_path)
    return run_shockcooling(toml, verbose)
end

end # module
