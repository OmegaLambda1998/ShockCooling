module Fitting

# External packages
using Unitful, UnitfulAstro
using KissMCMC 
using Distributions
using Supernovae

# Internal packages
using ..Models

# Exports
export run_mcmc

function get_prior(model::Model, numwalkers)
    ks = sort!(collect(keys(model.constraints)))
    x0 = [rand(model.constraints[k][1]) for (i, k) in enumerate(ks), j in 1:numwalkers]
    return x0
end

function prior_constraint(model::Model, param::Dict)
    rtn = 0
    for k in keys(model.constraints)
        model_min = minimum(model.constraints[k][1]) * model.constraints[k][2]
        model_max = maximum(model.constraints[k][1]) * model.constraints[k][2]
        if param[k] < model_min 
            return -Inf
        elseif param[k] > model_max
            return -Inf 
        end
        rtn += logpdf(model.constraints[k][1], ustrip(param[k]))
    end
    return rtn
end

function likelihood_function(model::Model, supernova::Supernova)
    function llhood(param)
        ks = sort!(collect(keys(model.constraints)))
        param = Dict(k => param[i] * model.constraints[k][2] for (i, k) in enumerate(ks)) 
        prior = prior_constraint(model, param)
        if prior == -Inf
            return prior
        end
        model_flux = run_model(model, param, supernova)
        chi = -0.5 * sum(((model_flux .- get_flux(supernova)) ./ get_flux_err(supernova)) .^ 2)
        return chi + prior
    end
    return llhood
end

function run_mcmc(config::Dict, model::Model, supernova::Supernova)
    @debug "Running with $(Threads.nthreads()) threads"
    numwalkers = config["numwalkers"] 
    thinning = config["thinning"]
    burnin = config["burnin"]
    numsamples_perwalker = config["numsamples_perwalker"]

    x0 = get_prior(model, numwalkers)
    llhood = likelihood_function(model, supernova)

    @info "Burning in"
    start_time = time()
    burn_chains, accept_ratio, llhoodvals, _ = emcee(llhood, x0; niter = burnin)
    @info "Finished in $((time() - start_time) / 60) mins"
    @info "Burn in has acceptance ratio of $accept_ratio"
    @info "Running"
    start_time = time()
    chains, accept_ratio, llhoodvals, _ = emcee(llhood, burn_chains[:, :, end]; niter = numsamples_perwalker)
    @info "Finished in $((time() - start_time) / 60) mins"
    @info "Run has acceptance ration of $accept_ratio"
    @info "MCMC Finished"

    return x0, burn_chains, chains, llhoodvals 
end

end
