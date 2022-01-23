module ShockCooling

# External packages
using ArgParse
using Supernovae
using TOML

# Internal Packages
include("Models.jl")
using .Models

# Exports
export run_mcmc

# Data function
function load_data(data_config::Dict)
    input_path = joinpath(dirname(@__DIR__), "Inputs")
    data_path = data_config["data_path"]
    if !isabspath(data_path)
        data_path = joinpath(input_path, data_path)
    end
    return Supernova(data_path)
end

# Model function
function load_model(model_config::Dict)
    constraint_dict = model_config["constraints"]
    constraints = get_constraints(constraint_dict)
    model = getfield(Models, Symbol(model_config["name"]))
    return model(constraints)
end

function load_model(model_list::Vector)
    return load_model.(model_list)
end

# Main function
function run_mcmc(toml::Dict)
    supernova = load_data(toml["data"])
    @show supernova.name
    models = load_model(toml["model"])
    @show models
end

function run_mcmc(toml_path::AbstractString)
    toml = TOML.parsefile(toml_path)
    return run_mcmc(toml)
end

function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--verbose", "-v"
            help = "Increase level of logging verbosity"
            action = :store_true
        "toml"
            help = "Path to toml input file"
            required = true
    end

    return parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    verbose = args["verbose"]
    toml_path = args["toml"]
    run_mcmc(toml_path)
end

end # module
