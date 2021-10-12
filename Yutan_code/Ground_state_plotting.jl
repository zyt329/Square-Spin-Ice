using JLD
using Plots

path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/"
name =["Square_Spin_Ice_Measurement_J2=1.0_L1=4_L2=6_hmin=0.01_hmax=10.0_Date_Sun_26_Sep_2021_03_42_26.jld",#=
"Square_Spin_Ice_Measurement_L=4_hmin=0.1_hmax=0.5_Date_Tue_27_Jul_2021_15_17_10.jld",
"Square_Spin_Ice_Measurement_tilted20sites_L1=11_L2=11_hmin=0.01_hmax=0.5_Date_Wed_04_Aug_2021_01_43_25.jld",
"Square_Spin_Ice_Measurement_J2=0.4444444444444444_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_01_46_36.jld",
"Square_Spin_Ice_Measurement_J2=0.6666666666666666_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_01_50_55.jld",
"Square_Spin_Ice_Measurement_J2=0.8888888888888888_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_01_55_58.jld",
"Square_Spin_Ice_Measurement_J2=1.1111111111111112_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_02_02_19.jld",
"Square_Spin_Ice_Measurement_J2=1.3333333333333333_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_02_07_59.jld",
"Square_Spin_Ice_Measurement_J2=1.5555555555555556_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_02_13_27.jld",
"Square_Spin_Ice_Measurement_J2=1.7777777777777777_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_02_18_46.jld",
"Square_Spin_Ice_Measurement_J2=2.0_L1=4_L2=4_hmin=0.01_hmax=0.5_Date_Thu_05_Aug_2021_02_24_00.jld"=#]


save_path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/plots/"

function plotting(savename::String, xlabel::String, Quantity::Array)
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
    plot(title = "m : J = 1, $(Quantity) sites lattice")
    for i in 1:length(name)
        h_vals = []
        m = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(m, result[3]...)

        plot!(
        h_vals,
        m,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "m",
        legend = :topright,
        label = "m ($(Quantity[1])=$(Quantity[i+1]))"
        )
    end
    savefig(save_path * "m_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")

    plot(title = "S_pi : J = 1, $(Quantity) sites lattice")
    for i in 1:length(name)
        h_vals = []
        S_pi = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(S_pi, result[4]...)

        plot!(
        h_vals,
        S_pi,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "S_pi",
        legend = :topright,
        label = "S_pi ($(Quantity[1])=$(Quantity[i+1]))"
        )

    end
    savefig(save_path * "S_pi_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")

    plot(title = "Fidelity : J = 1, $(Quantity) sites lattice")
    for i in 1:length(name)
        h_vals = []
        Fidelity = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(Fidelity, result[5]...)

        plot!(
        h_vals,
        Fidelity,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "Fidelity",
        legend = :topright,
        label = "Fidelity ($(Quantity[1])=$(Quantity[i+1]))"
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
    plot(title = "S_entangle : J = 1, $(Quantity) sites lattice")
    for i in 1:length(name)
        h_vals = []
        S_entangle_vals = []
        result = load(path*name[i], "result")
        push!(h_vals, result[1]...)
        push!(S_entangle_vals, result[6]...)

        plot!(
        h_vals,
        S_entangle_vals,
        dpi = 800,
        xlabel = xlabel,
        ylabel = "S_entangle",
        legend = :topright,
        label = "S_entangle ($(Quantity[1])=$(Quantity[i+1]))"
        )
    end
    savefig(save_path * "S_entangle_" *savename*"hmin=$(h_vals[1])_hmax=$(h_vals[end])"* ".png")
end
plotting("Square_Spin_Ice_Measurement_24sites_Spi_0_Spi_compare", "h", vcat(["N"],[24]))
#println(Fidelity)
