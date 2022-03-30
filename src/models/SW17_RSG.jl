# Based of Sapir & Waxman 2017
# DOI: 10.3847/1538-4357/aa64df
export SW17_RSG
mutable struct SW17_RSG <: Model
    name::AbstractString
    class::AbstractString
    parameter_names::Dict{String, LaTeXString}
    constraints::Dict
end

# Default parameter names
function SW17_RSG(name::AbstractString, class::AbstractString, constraints::Dict)
    parameter_names = Dict("R" => L"R_{e}~[R_{\odot}]", "M" => L"M_{e}~[M_{\odot}]", "v" => L"v_{s}~[\frac{km}{s}]", "t" => L"t_{off}~[Days]")
    return SW17_RSG(name, class, parameter_names, constraints)
end

# Default name and parameter names
function SW17_RSG(constraints::Dict)
    name = "Sapir & Waxman (2017) RSG"
    class = "SW17_RSG"
    return SW17_RSG(name, class, constraints)
end

function bolometric_luminosity(model::SW17_RSG, param::Dict, observation::Observation)
    # Load parameters
    r = param["R"] # Envelope radius
    m = param["M"] # Envelope mass
    v = param["v"] # Envelope velocity
    t = (observation.time - param["t"]) # Time since explosion
    t = max(0 * u"d", t)
    if t == 0 * u"d"
        return 0 * u"erg / s"
    end

    # Physical values
    k = 1 # 0.34 cm^2 / g
    mc = 1u"Msun"

    # Rescaled input parameters
    menv = m + mc
    vs = v / (10^8.5 * u"cm / s") |> NoUnits
    M0 = m / 1u"Msun" |> NoUnits
    td = t / 1u"d" |> NoUnits
    R13 = r / (10^13 * u"cm") |> NoUnits

    # n dependent values
    α = 0.8
    a = 1.67
    ϵ2 = -0.086
    fp = (m / mc) ^ 0.5 
    Lu = (1.88e42) * u"erg / s"

    # Other values
    ttr = 19.5u"d" * (k * M0 / vs) ^ 0.5

    LRW = Lu * vs * vs * R13 * ((vs * td * td / (fp * (menv / 1u"Msun") * k)) ^ ϵ2) / k
    e = -(a * t / ttr) ^ α

    L = LRW * exp(e)

    return L
end

function radius(model::SW17_RSG, param::Dict, observation::Observation)
    L = bolometric_luminosity(model, param, observation)
    T = temperature(model, param, observation)
    σ = 5.6704e-8 * u"W / m^2 / K^4" # Stefan Boltzmann constant

    R = (L / (4π * σ * T ^ 4)) ^ 0.5

    return R
end

function temperature(model::SW17_RSG, param::Dict, observation::Observation)
    r = param["R"] # Envelope radius
    v = param["v"] # Envelope velocity
    m = param["M"] # Envelope mass
    t = observation.time - param["t"] # Time since explosion
    t = max(0 * u"d", t)

    # Physical values
    k = 1
    mc = 1u"Msun"
    menv = m + mc

    # Rescaled input parameters
    vs = v / (10^8.5 * u"cm / s") |> NoUnits
    td = t / 1u"d" |> NoUnits

    M0 = m / 1u"Msun" |> NoUnits

    R13 = r / (10^13 * u"cm") |> NoUnits


    # n dependent values
    Tu = 2.05e4 * u"K"
    fp = (m / mc) ^ 0.5 
    ϵ1 = 0.027

    T = Tu * ((vs * vs * td * td / (fp * (menv / 1u"Msun") * k)) ^ ϵ1) * ((R13 / (td * td * k)) ^ 0.25)
    return T
end

function model_flux(model::SW17_RSG, param::Dict, observation::Observation)
    T = temperature(model, param, observation)
    return synthetic_flux(observation.filter, T)
end

function run_model(model::SW17_RSG, param::Dict, supernova::Supernova)
    m_flux = [model_flux(model, param, obs) for obs in supernova.lightcurve.observations] .|> u"erg / s / cm^2 / Hz"
    dist = 10u"pc"
    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    abs_mag = @. -48.6 - 2.5 * (log10(ustrip(m_flux)) + log10(uconvert(NoUnits, R / dist)^2))
    replace!(abs_mag, NaN => -10)
    return abs_mag * u"AB_mag"
end
