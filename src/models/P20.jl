# Based of Piro 2020
# DOI: 10.3847/1538-4357/abe2b1

export P20
mutable struct P20 <: Model
    name::AbstractString
    parameter_names::Dict{String, LaTeXString}
    constraints::Dict
end

# Default parameter names
function P20(name::AbstractString, constraints::Dict)
    parameter_names = Dict("R" => L"R_{e}~[R_{\odot}]", "M" => L"M_{e}~[M_{\odot}]", "v" => L"v_{e}~[\frac{km}{s}]", "t" => L"t_{off}~[Days]")
    return P20(name, parameter_names, constraints)
end

# Default name and parameter names
function P20(constraints::Dict)
    name = "Piro (2020)"
    return P20(name, constraints)
end

function bolometric_luminosity(model::P20, param::Dict, observation::Observation)
    #TODO fix up the equation to not rely on weird numerical values
    # Load parameters
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["v"] # Envelope velocity
    t = (observation.time - param["t"]) # Time since explosion
    t = max(1e-10u"d", t)
    
    # All other values
    k = 0.34 * u"cm^2 / g"
    c = 299792458 * u"m / s" # Speed of light
    n = 10
    δ = 1.1
    K = (n - 3) * (3 - δ) / (4π * (n - δ))

    td = (3 * k * K * m / ((n - 1) * v * c)) ^ 0.5

    L0 = π * (n - 1) * c * r * v * v / (3 * (n - 5) * k) 

    if t < td
        L = L0 * (td / t) ^ (4 / (n - 2))
    else
        e = -0.5 * ((t * t / (td * td)) - 1)
        L = L0 * exp(e)
    end
    return L
end

function radius(model::P20, param::Dict, observation::Observation)
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["v"] # Envelope velocity
    t = observation.time - param["t"] # Time since explosion
    t = max(1e-10u"d", t)

    # All other values
    k = 0.34 * u"cm^2 / g"
    n = 10
    δ = 1.1
    K = (n - 3) * (3 - δ) / (4π * (n - δ))


    tph = (3 * k * K * m / (2 * v * v * (n - 1))) ^ 0.5

    if t < tph
        R = t * v * (tph / t) ^ (2 / (n - 1))
    else
        R = t * v * (1 + ((δ - 1) / (n - 1)) * ((t * t / (tph * tph)) - 1)) ^ (-1 / (δ - 1))
    end

    return R
end

function temperature(model::P20, param::Dict, observation::Observation)
    L = bolometric_luminosity(model, param, observation) |> u"erg / s"
    R = radius(model, param, observation) |> u"km"
    σ = 5.6704e-8 * u"W / m^2 / K^4" # Stefan Boltzmann constant
    T = (L / (4π * R * R * σ)) ^ 0.25 |> u"K"
    return T
end

function model_flux(model::P20, param::Dict, observation::Observation)
    T = temperature(model, param, observation)
    return synthetic_flux(observation.filter, T)
end

function run_model(model::P20, param::Dict, supernova::Supernova)
    m_flux = [model_flux(model, param, obs) for obs in supernova.lightcurve.observations] .|> u"erg / s / cm^2 / Hz"
    dist = 10u"pc"
    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    abs_mag = @. -48.6 - 2.5 * (log10(ustrip(m_flux)) + log10(uconvert(NoUnits, R / dist)^2))
    return abs_mag * u"AB_mag"
end
