module Models

# External Packages
using Supernovae
using LaTeXStrings
using Unitful, UnitfulAstro
using Distributions

# Exports
export Model
export get_model
export run_model
export get_constraints
export bolometric_luminosity, radius, temperature

abstract type Model end

# Load in all models
# All models must export themselves
model_path = joinpath(@__DIR__, "models")
for path in readdir(model_path)
    if isfile(joinpath(model_path, path))
        include(joinpath(model_path, path))
    end
end

function get_model(model_name::String)
    return getfield(Models, Symbol(model_name))
end

function run_model(model::Model, param::Dict, observation::Observation)
    error("generic run_model is being used, no method specified for model $(typeof(model))")
end

function get_constraints(constraints_dict::Dict)
    constraints = Dict()
    for param in keys(constraints_dict)
        param_dict = constraints_dict[param]
        param_unit = uparse(param_dict["unit"], unit_context = [Unitful, UnitfulAstro])
        param_prior = getfield(Distributions, Symbol(param_dict["prior"]))
        param_values = param_dict["values"]
        param_min = get(param_dict, "min", -Inf)
        param_max = get(param_dict, "max", Inf)
        if param_prior == Normal
            prior = TruncatedNormal(param_values..., param_min, param_max)
        else
            prior = Truncated(param_prior(param_values...), param_min, param_max)
        end
        constraints[param] = (prior, param_unit)
    end
    return constraints
end



end
