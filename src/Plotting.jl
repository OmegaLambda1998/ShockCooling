module Plotting

# External Packages
using CairoMakie
CairoMakie.activate!(type = "svg")
using Unitful, UnitfulAstro
using Supernovae
using Random
Random.seed!(0)
using KernelDensity
using LaTeXStrings
using Colors


# Internal Packages
using ..Models
using ..Fitting

# Exports
export plot_luminosity, plot_luminosity!
export plot_radius, plot_radius!
export plot_temperature, plot_temperature!
export plot_contour, plot_contour!
export plot_comparison, plot_comparison!

# Markers
markers_labels = shuffle([
    :rect,
    :star5,
    :diamond,
    :hexagon,
    :cross,
    :xcross,
    :utriangle,
    :dtriangle,
    :ltriangle,
    :rtriangle,
    :pentagon,
    :star4,
    :star8,
    :vline,
    :hline,
    :x,
    :+,
    :circle
   ])

# Colours
colour_labels = shuffle(["salmon", "coral", "tomato", "firebrick", "crimson", "red", "orange", "green", "forestgreen", "seagreen", "olive", "lime", "charteuse", "teal", "turquoise", "cyan", "navyblue", "midnightblue", "indigo", "royalblue", "slateblue", "steelblue", "blue", "purple", "orchid", "magenta", "maroon", "hotpink", "deeppink", "saddlebrown", "brown", "peru", "tan"])

function plot_temperature(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    temp_unit = uparse(get(units, "temp", "K"), unit_context = [Unitful, UnitfulAstro])
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
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    temp_unit = uparse(get(units, "temp", "K"), unit_context = [Unitful, UnitfulAstro])

    T = [temperature(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = ustrip.(T .|> temp_unit) .+ 1e-10
    scatter!(ax, x, y, label=model.name)
end

function plot_radius(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    rad_unit= uparse(get(units, "radius", "Rsun"), unit_context = [Unitful, UnitfulAstro])

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
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    rad_unit = uparse(get(units, "radius", "Rsun"), unit_context = [Unitful, UnitfulAstro])

    R = [radius(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = ustrip.(R .|> rad_unit)
    scatter!(ax, x, y, label=model.name)
end

function plot_luminosity!(fig, ax, model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    lum_unit = uparse(get(units, "luminosity", "erg / s"), unit_context = [Unitful, UnitfulAstro])

    L = [bolometric_luminosity(model, param, obs) for obs in supernova.lightcurve.observations]
    x = ustrip.((get(supernova, "time") .- param["t"]) .|> time_unit)
    y = log10.(ustrip.(L .|> lum_unit))
    scatter!(ax, x, y, label=model.name)
end

function plot_luminosity(model::Model, param::Dict, supernova::Supernova, plot_config::Dict)
    fig = Figure()
    units = get(plot_config, "unit", Dict())
    time_unit = get(units, "time", nothing)
    if isnothing(time_unit)
        time_unit = unit(get(supernova, "time")[1])
    else
        time_unit = uparse(time_unit, unit_context = [Unitful, UnitfulAstro])
    end
    lum_unit = uparse(get(units, "luminosity", "erg / s"), unit_context = [Unitful, UnitfulAstro])

    path = get(plot_config, "path", nothing)
    ax = Axis(fig[1, 1], xlabel = "Time [$time_unit]", ylabel = "Luminosity [$lum_unit]", title = "Bolometric Luminosity") 
    plot_luminosity!(fig, ax, model, param, supernova, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, ax
end

function plot_contour!(fig, gax, model::Model, chain, plot_config::Dict)
    ks = sort!(collect(keys(model.constraints)))
    # Tranpose chain so each element is a parameter rather than a walker
    chain = collect(eachrow(reduce(hcat, chain)))
    n = length(ks)
    @debug "Creating a $(n)x$(n) sized contour plot"
    for i in 1:n
        for j in 1:i
            ax = Axis(gax[i, j], xlabel = model.parameter_names[ks[j]], ylabel = model.parameter_names[ks[i]])
            if i == j
                density!(ax, chain[i], color = "darkblue")
                hist!(ax, chain[i], color = :transparent, normalization = :pdf, strokewidth = 1, strokecolor = "black")
                ylims!(ax, low = 0)
                xlims!(ax, minimum(chain[i]), maximum(chain[i]))
                hideydecorations!(ax)
                if i != n
                    hidexdecorations!(ax)
                end
            else
                k = kde(hcat(chain[j], chain[i]))
                contour!(ax, k.x, k.y, k.density)
                if j != 1
                    hideydecorations!(ax, grid=false)
                end
                if i != n
                    hidexdecorations!(ax, grid=false)
                end
                xlims!(ax, minimum(chain[j]), maximum(chain[j]))
                ylims!(ax, minimum(chain[i]), maximum(chain[i]))
            end
        end
    end
end

function plot_contour(model::Model, chain, plot_config::Dict)
    fig = Figure(resolution=(1600, 1200))
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_contour!(fig, gax, model, chain, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

function plot_comparison!(fig, gax, model::Model, supernova::Supernova, param::Dict, chain, plot_config::Dict)
    data_type = get(plot_config, "data_type", "flux")
    @debug "Plotting data type set to $data_type"
    units = get(plot_config, "unit", Dict())
    time_unit = uparse(get(units, "time", "d"), unit_context = [Unitful, UnitfulAstro])
    names = get(plot_config, "names", nothing)
    @debug "Generating all plot vectors"
    filters = Set([obs.filter.name for obs in supernova.lightcurve.observations])
    lc_ax = Axis(gax[1, 1], xlabel = "Time [$time_unit]")
    res_ax = Axis(gax[2, 1], xlabel = "Time [$time_unit]")
    if data_type in ["magnitude", "abs_magnitude"]
        lc_ax.yreversed = true
        res_ax.yreversed = true
    end
    colours, markers = plot_lightcurve!(fig, lc_ax, supernova, plot_config)
    num_chains = min(1000, length(chain))
    chains = []
    for c in shuffle(chain)[1:num_chains]
        d = Dict()
        for (i, k) in enumerate(sort!(collect(keys(model.parameter_names))))
            d[k] = c[i] * model.constraints[k][2]
        end 
        push!(chains, d)
    end
    for filt in filters
        sn = filter(obs -> obs.filter.name == filt, supernova)
        time = get(sn, "time")
        m_absmag = run_model(model, param, sn)
        c_absmag = [run_model(model, c, sn) for c in chains]
        m_mag = absmag_to_mag.(m_absmag, sn.redshift)
        c_mag = [absmag_to_mag.(c, sn.redshift) for c in c_absmag]
        m_flux = mag_to_flux.(m_mag, sn.zeropoint)
        c_flux = [mag_to_flux.(c, sn.zeropoint) for c in c_mag]
        if data_type == "flux"
            data_unit = uparse(get(units, "data", "µJy"), unit_context = [Unitful, UnitfulAstro])
            m_data = m_flux
            c_data = c_flux
        elseif data_type == "magnitude"
            data_unit = uparse(get(units, "data", "AB_mag"), unit_context = [Unitful, UnitfulAstro])
            m_data = m_mag
            c_data = c_mag
        elseif data_type == "abs_magnitude"
            data_unit = uparse(get(units, "data", "AB_mag"), unit_context = [Unitful, UnitfulAstro])
            m_data = m_absmag
            c_data = c_absmag
            lc_ax.yreversed = true
            res_ax.yreversed = true
        else
            error("Unknown data type: $data_type. Possible data types are [flux, magnitude, abs_magnitude]")
        end
        lc_ax.ylabel = "$data_type [$data_unit]"
        res_ax.ylabel = "model - data [$data_unit]"
        data = get(sn, data_type)
        data_err = get(sn, "$(data_type)_err")
        lines!(lc_ax, ustrip(time), ustrip(m_data .|> data_unit), color = colours[filt])
        for c in c_data
            lines!(lc_ax, ustrip(time), ustrip(c .|> data_unit), color = alphacolor(parse(Colorant, colours[filt]), 0.01))
        end
        lines!(res_ax, ustrip(time), [0 for t in time], color = "black")
        scatter!(res_ax, ustrip(time), ustrip((m_data .|> data_unit)) .- ustrip((data .|> data_unit)), color = colours[filt])
        errorbars!(res_ax, ustrip(time), ustrip((m_data .|> data_unit)) .- ustrip((data .|> data_unit)), ustrip(data_err .|> data_unit), color = colours[filt])
    end
end

function plot_comparison(model::Model, supernova::Supernova, param::Dict, chain, plot_config::Dict)
    fig = Figure()
    path = get(plot_config, "path", nothing)
    gax = fig[1, 1] = GridLayout()
    plot_comparison!(fig, gax, model, supernova, param, chain, plot_config)
    if !isnothing(path)
        save(path, fig)
    end
    return fig, gax
end

end
