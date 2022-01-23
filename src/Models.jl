module Models

# External Packages
using Supernovae
using LaTeXStrings
using Unitful, UnitfulAstro
using Distributions

# Exports
export run_model
export get_constraints

abstract type Model end

# Load in all models
# Models should export themselves
model_path = joinpath(@__DIR__, "models")
for path in readdir(model_path)
    if isfile(joinpath(model_path, path))
        include(joinpath(model_path, path))
    end
end

function run_model(model::Model, param::Dict, observation::Observation)
    error("generic run_model is being used, no method specified for model $(typeof(model))")
end

function get_constraints(constraint_dict::Dict)
    constraints = Dict()
    for param in keys(constraints_dict)
        param_dict = constraints_dict[param]
        param_unit = uparse(param_dict["unit"])
        param_prior = getfield(Distributions, Symbol(param_dict["prior"]))
        param_values = param_dict["values"] .* param_unit
        param_min = get(param_dict, "min", -Inf) .* param_unit
        param_max = get(param_dict, "max", Inf) .* param_unit
        if param_prior == Normal
            prior = TruncatedNormal(param_values..., param_min, param_max)
        else
            prior = TruncatedNormal(param_prior(param_values...), param_min, param_max)
        end
        constraints[param] = prior
    end
    return constraints
end

end
