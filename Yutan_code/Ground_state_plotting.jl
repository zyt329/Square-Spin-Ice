using JLD
using Plots

path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/"
name =["Square_Spin_Ice_Measurement_tilted20sites_L1=11_L2=11_hmin=0.01_hmax=0.5_Date_Wed_04_Aug_2021_01_43_25.jld","Square_Spin_Ice_Measurement_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Sun_01_Aug_2021_01_33_21.jld"]


save_path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/plots/"

function plotting(savename::String, xlabel::String, N::Array)
    h_vals = []
    m = []
    S_pi = []
    Fidelity = []
    #Fcl_vals = []
    #Fqm_vals = []
    S_entangle_vals = []
    #push!(S_pi, result[3]...)
    #push!(Fidelity, result[4]...)
    #push!(Fcl_vals, result[5]...)
    #push!(Fqm_vals, result[6]...)
    #push!(S_entangle_vals, result[7]...)
    plot()
    for i in 1:length(name)
        h_vals = []
        m = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(m, result[2]...)

        plot!(
        h_vals,
        m,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "m",
        title = "m : J = 1, $(N) sites lattice",
        legend = :topright,
        label = "m ($(N[i]))"
        )
    end
    savefig(save_path * "m_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")

    plot()
    for i in 1:length(name)
        h_vals = []
        S_pi = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(S_pi, result[3]...)

        plot!(
        h_vals,
        S_pi,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "S_pi",
        title = "S_pi : J = 1, $(N) sites lattice",
        legend = :topright,
        label = "S_pi ($(N[i]))"
        )

    end
    savefig(save_path * "S_pi_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")

    plot()
    for i in 1:length(name)
        h_vals = []
        Fidelity = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(Fidelity, result[4]...)

        plot!(
        h_vals,
        Fidelity,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "Fidelity",
        title = "Fidelity : J = 1, $(N) sites lattice",
        legend = :topright,
        label = "Fidelity ($(N[i]))"
        )

    end
    savefig(save_path * "Fidelity_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")
#=
    plot(
    h_vals,
    Fcl_vals,
    dpi = 800,
    xlabel = xlabel,
    ylabel = "Fcl",
    title = "Fcl : J = 1, $(N) sites lattice",
    legend = :topright,
    label = "Fcl"
    )
    savefig(save_path * "Fcl_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")

    plot(
    h_vals,
    Fqm_vals,
    dpi = 800,
    xlabel = xlabel,
    ylabel = "Fqm",
    title = "Fqm : J = 1, $(N) sites lattice",
    legend = :topright,
    label = "Fqm"
    )
    savefig(save_path * "Fqm_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")=#
    plot()
    for i in 1:length(name)
        h_vals = []
        S_entangle_vals = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        if i == 1
            push!(S_entangle_vals, result[5]...)
        else
            push!(S_entangle_vals, result[7]...)
        end

        plot!(
        h_vals,
        S_entangle_vals,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "S_entangle",
        title = "S_entangle : J = 1, $(N) sites lattice",
        legend = :topright,
        label = "S_entangle ($(N[i]))"
        )
    end
    savefig(save_path * "S_entangle_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")
end
plotting("Square_Spin_Ice_Measurement_20sites", "h", [20, 16])
#println(Fidelity)
