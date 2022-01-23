# Based of Piro 2015
# DOI: 10.1088/2041-8205/808/2/L51 
#
# Using the form found in Arcavi Et. Al 2017
# DOI: 10.3847/2041-8213/aa5be1

export P15
mutable struct P15 <: Model
    parameter_names::Dict{String, LaTeXString}
    constraints::Dict
end

# Default parameter names
function P15(constraints::Dict)
    parameter_names = Dict("R" => L"R_{e}~[R_{\odot}]", "M" => L"M_{e}~[M_{\odot}]", "v" => L"v_{e}~[\frac{km}{s}]", "t" => L"t_{off}~[Days]")
    return P15(parameter_names, constraints)
end

function bolometric_luminosity(model::P15, param::Dict, observation::Observation)
    #TODO fix up the equation to not rely on weird numerical values
    # Load parameters
    r = 1e-13 * param["R"] |> u"cm" # Envelope radius in 10^13 cm
    m = u"Msun" |> param["M"] # Envelope mass in solar mass
    v = 1e-9 * param["v"] |> u"cm / s" # Envelope velocity in 10^9 cm / s
    t = (obvservation.time - param["t"]) |> u"s" # Time since explosion in s
    
    # All other values
    k = 0.34 * u"cm^2 / g" # Opacity
    mc = 1 # Core mass. Only ever present as mc / Msun so no need to specify units

    # Numerical factors
    c1 = 8.27e42 * u"cm"
    c2 = -4.134e-11 * u"cm / g / s"
    c3 = 2e4 # Unitless
    c4 = 0.01 * u"Msun"

    # Exponential Factor
    exp = (c2 / k) * t * ((t * v) + (c3 * r)) * (mc ^ 0.01) * (c4 / m)

    # L at t=0
    L0 = (c1 / k) * v * v * r * (mc ^ 0.01)

    # Bolometric Luminosity
    L = L0 * exp

    return L
end

function radius(model::P15, param::Dict, observation::Observation)
    r = param["R"] # Envelope radius
    v = param["v"] # Envelope velocity
    t = observation.time - param["t"] # Time since explosion

    R = r + v * t
    
    return R
end

function temperature(model::P15, param::Dict, observation::Observation)
    L = bolometric_luminosity(model, param, observation)
    R = radius(model, param, observation)
    σ = 5.67e-8 * u"W / m^2 / K^4" # Stefan Boltzmann constant
    T = (L / (4π * R * R * σ)) ^ 0.25
end

function flux(model::P15, param::Dict, observation::Observation)
    T = temperature(model, param, observation)
    return flux(observation.filter, T)
end

function run_model(model::P15, param::Dict, supernova::Supernova)
    model_flux = flux.(model, param, supernova.lightcurve.observations)
    flux_unit = unit(supernovae.lightcurve.observations[1].flux)
    return model_flux .|> flux_unit
end
