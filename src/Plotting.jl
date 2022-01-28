module Plotting

# External Packages
using CairoMakie
CairoMakie.activate!(type = "svg")
using Unitful, UnitfulAstro
using Supernovae


# Internal Packages
using ..Models
using ..Fitting

# Exports
export plot_luminosity, plot_luminosity!
export plot_radius, plot_radius!
export plot_temperature, plot_temperature!
export plot_prior, plot_prior!
export plot_walkers, plot_walkers!
export plot_contour, plot_contour!
export plot_comparison, plot_comparison!

function plot_temperature(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    temp_unit = uparse(get(units, "temp", "K"))
    path = get(plot_config, "path", nothing)
    ax = Axis(fig[1, 1], yscale = log10, xlabel = "Time [$time_unit]", ylabel = "Temperature [$temp_unit]", title = "Temperature plot") 
    plot_temperature!(fig, ax, model, param, supernova, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, ax
end

function plot_temperature!(fig, ax, model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    temp_unit = uparse(get(units, "temp", "K"))

    T = [temperature(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = ustrip.(T .|> temp_unit)
    lines!(ax, x, y, label=model.name)
end

function plot_radius(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    rad_unit= uparse(get(units, "radius", "Rsun"))

    path = get(plot_config, "path", nothing)
    ax = Axis(fig[1, 1], yscale = log10, xlabel = "Time [$time_unit]", ylabel = "Radius [$rad_unit]", title = "Radius") 
    plot_radius!(fig, ax, model, param, supernova, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, ax
end

function plot_radius!(fig, ax, model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    rad_unit = uparse(get(units, "radius", "Rsun"))

    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = ustrip.(R .|> rad_unit)
    lines!(ax, x, y, label=model.name)
end

function plot_luminosity!(fig, ax, model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    lum_unit = uparse(get(units, "luminosity", "erg / s"))

    L = [bolometric_luminosity(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = log10.(ustrip.(L .|> lum_unit))
    lines!(ax, x, y, label=model.name)
end

function plot_luminosity(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit)
    end
    lum_unit = uparse(get(units, "luminosity", "erg / s"))

    path = get(plot_config, "path", nothing)
    ax = Axis(fig[1, 1], xlabel = "Time [$time_unit]", ylabel = "Luminosity [$lum_unit]", title = "Bolometric Luminosity") 
    plot_luminosity!(fig, ax, model, param, supernova, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, ax
end

function plot_prior!(fig, gax, model::Model, priors, plot_config::Dict)
    ks = sort!(collect(keys(model.constraints)))
    for (i, k) in enumerate(ks)
        ax = Axis(gax[i, 1], xlabel = "$(model.parameter_names[k])")
        hist!(ax, priors[i, :])
    end
end

function plot_prior(model::Model, priors, plot_config::Dict)
    fig = Figure()
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_prior!(fig, gax, model, priors, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

function plot_walkers!(fig, gax, model::Model, walkers, plot_config::Dict)
    ks = sort!(collect(keys(model.constraints)))
    dim, nwalker, steps = size(walkers)
    for (i, k) in enumerate(ks)
        ax = Axis(gax[i, 1], ylabel = model.parameter_names[k], xlabel="Steps")
        for j in 1:nwalker
            chain = walkers[i, j, :]
            lines!(ax, 1:steps, chain, color=:grey)
        end
    end
end

function plot_walkers(model::Model, walkers, plot_config::Dict)
    fig = Figure()
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_walkers!(fig, gax, model, walkers, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

function plot_contour!(fig, gax, model::Model, chain, llhood, plot_config::Dict)
    ks = sort!(collect(keys(model.constraints)))
    dim, nwalker, steps = size(chain)
    llhood = vec(reshape(llhood, (1, nwalker * steps)))
    zs = zs * zs'
    for (i, ki) in enumerate(ks)
        for (j, kj) in enumerate(ks)
            if j > i
                continue
            end
            ax = Axis(fig[i, j])
            if i == 1
                ax.ylabel = model.parameter_names[ki]
            end
            if j == dim
                ax.xlabel = model.parameter_names[kj]
            end
            xs = vec(reshape(chain[i, :, :], (1, nwalker * steps)))
            if i == j
                density!(ax, xs)
            else
                ys = vec(reshape(chain[j, :, :], (1, nwalker * steps)))
                contour!(ax, xs, ys, zs)
            end
        end
    end
end

function plot_contour(model::Model, chain, llhood, plot_config::Dict)
    fig = Figure()
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_contour!(fig, gax, model, chain, llhood, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

function plot_comparison!(fig, gax, model::Model, supernova::Supernova, param::Dict, plot_config::Dict)
    time = get(supernova, "time")
    data_type = get(plot_config, "data_type", "flux")
    if data_type == "flux"
        data = get(supernova, "flux")
    elseif data_type == "magnitude"
        data = get(supernova, "magnitude")
    elseif data_type == "abs_magnitude"
        data = get(supernova, "abs_magnitude")
    else
        error("Unknown data type: $data_type. Possible data types are [flux, magnitude, abs_magnitude]")
    end
    m_absmag = run_model(model, param, supernova)
    m_mag = absmag_to_mag.(m_absmag, supernova.redshift)
    m_flux = mag_to_flux.(m_mag, supernova.zeropoint)
    lc_ax = Axis(gax[1, 1])
    res_ax = Axis(gax[2, 1])
    if data_type == "flux"
        @info unit(data[1])
        m_data = m_flux .|> u"µJy" 
    elseif data_type == "magnitude"
        m_data = m_mag .|> u"AB_mag" 
        lc_ax.yreversed = true
        res_ax.yreversed = true
    elseif data_type == "abs_magnitude"
        m_data = m_absmag .|> u"AB_mag" 
        lc_ax.yreversed = true
        res_ax.yreversed = true
    else
        error("Unknown data type: $data_type. Possible data types are [flux, magnitude, abs_magnitude]")
    end
    plot_lightcurve!(fig, lc_ax, supernova, plot_config)
    scatter!(lc_ax, ustrip(time), ustrip(m_data))
    lines!(res_ax, ustrip(time), [0 for t in time])
    scatter!(res_ax, ustrip(time), ustrip(m_data) .- ustrip(data))
end

function plot_comparison(model::Model, supernova::Supernova, param::Dict, plot_config::Dict)
    fig = Figure()
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_comparison!(fig, gax, model, supernova, param, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

end
