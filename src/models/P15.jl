# Based of Piro 2015
# DOI: 10.1088/2041-8205/808/2/L51 

export P15
mutable struct P15 <: Model
    name::AbstractString
    parameter_names::Dict{String, LaTeXString}
    constraints::Dict
end

# Default parameter names
function P15(name::AbstractString, constraints::Dict)
    parameter_names = Dict("R" => L"R_{e}~[R_{\odot}]", "M" => L"M_{e}~[M_{\odot}]", "v" => L"v_{e}~[\frac{km}{s}]", "t" => L"t_{off}~[Days]")
    return P15(name, parameter_names, constraints)
end

# Default name and parameter names
function P15(constraints::Dict)
    name = "Piro (2015)"
    return P15(name, constraints)
end

function bolometric_luminosity(model::P15, param::Dict, observation::Observation)
    # Load parameters
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["v"] # Envelope velocity
    t = (observation.time - param["t"]) # Time since explosion
    t = max(0 * u"d", t)
    
    # All other values
    k = 1 # Opacity in units of 0.34 cm^2 / g. Only every present as k / 0.34 cm^2/g so no need to specify units
    mc = 1 # Core mass. Only ever present as mc / Msun so no need to specify units
    c = 299792458 * u"m / s" # Speed of light
    mu = 0.01 * u"Msun" # Mass scaling
    vu = 2e9 * u"cm / s" # Velocity scaling
    Eu = 4e49 * u"erg" # Energy scaling
    tu = 0.9 * u"d" # Time scaling

    te = r / v
    E51 = ((v / vu) * (mc ^ 0.35) * ((m / mu) ^ 0.15)) ^ 2
    Ee = Eu * E51 * (mc ^ -0.7)  * ((m / mu) ^ 0.7)
    tp = tu * (k ^ 0.5) * (E51 ^ -0.25) * (mc ^ 0.17) * ((m / mu) ^ 0.57)

    L0 = te * Ee / (tp * tp)
    e = - t * (t + 2 * te) / (2 * tp * tp)

    L = L0 * exp(e)
    return L
end

function radius(model::P15, param::Dict, observation::Observation)
    r = param["R"] # Envelope radius
    v = param["v"] # Envelope velocity
    t = observation.time - param["t"] # Time since explosion
    t = max(0 * u"d", t)

    R = r + v * t
    return R
end

function temperature(model::P15, param::Dict, observation::Observation)
    L = bolometric_luminosity(model, param, observation)
    R = radius(model, param, observation)
    σ = 5.6704e-8 * u"W / m^2 / K^4" # Stefan Boltzmann constant
    T = (L / (4π * R * R * σ)) ^ 0.25
    return T
end

function model_flux(model::P15, param::Dict, observation::Observation)
    T = temperature(model, param, observation)
    return synthetic_flux(observation.filter, T)
end

function run_model(model::P15, param::Dict, supernova::Supernova)
    m_flux = [model_flux(model, param, obs) for obs in supernova.lightcurve.observations] .|> u"erg / s / cm^2 / Hz"
    dist = 10u"pc"
    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    abs_mag = @. -48.6 - 2.5 * (log10(ustrip(m_flux)) + log10(uconvert(NoUnits, R / dist)^2))
    return abs_mag * u"AB_mag"
end
