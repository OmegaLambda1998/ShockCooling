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
export save_chain, load_chain
export likelihood_function

function get_prior(model::Model, numwalkers)
    ks = sort!(collect(keys(model.constraints)))
    x0 = [[rand(model.constraints[k][1]) for k in ks] for j in 1:numwalkers]
    return x0
end

function prior_constraint(model::Model, param::Dict)
    return sum([logpdf(model.constraints[k][1], ustrip(param[k])) for k in keys(model.constraints)])
end

function likelihood_function(model::Model, supernova::Supernova)
    function llhood(param)
        ks = sort!(collect(keys(model.constraints)))
        param = Dict(k => param[i] * model.constraints[k][2] for (i, k) in enumerate(ks)) 
        prior = prior_constraint(model, param)
        # If outside physical bounds, don't bother evaluating the probability
        if !isfinite(prior)
            return prior
        end
        m_absmag = run_model(model, param, supernova)
        m_mag = absmag_to_mag.(m_absmag, supernova.redshift)
        m_flux = mag_to_flux.(m_mag, supernova.zeropoint) 
        chi = sum(((m_flux .- get(supernova, "flux")) ./ get(supernova, "flux_err")) .^ 2)
        r = length(m_flux) - length(ks)
        r_chi = -0.5 * chi / r
        return prior + r_chi 
    end
    return llhood
end

function save_chain(path::AbstractString, chain)
    chain = join([join(c, ",") for c in chain], "\n")
    open(path, "w") do io
        write(io, chain)
    end
end

function load_chain(path::AbstractString)
    chain = open(path, "r") do io
        r = readlines(io)
        chain = [[parse(Float64, s) for s in split(line, ",")] for line in r]
    end
    return chain
end

function run_mcmc(config::Dict, model::Model, supernova::Supernova)
    @debug "Running with $(Threads.nthreads()) threads"
    # Important config
    numwalkers = config["numwalkers"] 
    @debug "Numwalkers: $numwalkers"
    thinning = get(config, "thinning", 1)
    burnin = config["burnin"]
    @debug "Burnin: $burnin"
    iterations = config["iterations"]
    @debug "Iterations: $iterations"
    llhood = likelihood_function(model, supernova)

    # Unimportant config
    use_progress_meter = get(config, "progress_meter", true)

    @info "Generating priors"
    x0 = get_prior(model, numwalkers)

    @info "Running MCMC"
    chain, accept_ratio, logdensities, blob = emcee(llhood, x0; niter = iterations, nburnin = burnin, nthin = thinning, use_progress_meter = use_progress_meter)
    @info "MCMC Finished"

    chain, accept_ratio, logdensities, blob = squash_walkers(chain, accept_ratio, logdensities, blob)
    @info "MCMC has accept ratio of $accept_ratio"

    return x0, chain, accept_ratio, logdensities, blob 
end

end
