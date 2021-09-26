using JLD
using Plots

save_path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/plots/"

files = ["E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_J2=1_L1=4_L2=6_hmin=0.01_hmax=0.5_Date_Tue_10_Aug_2021_18_09_52.jld",
"E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_tilted20sites_L1=11_L2=11_hmin=0.01_hmax=0.5_Date_Tue_17_Aug_2021_19_44_55.jld",
"E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_J2=1_L1=4_L2=4_hmin=0.01_hmax=1.0_Date_Wed_18_Aug_2021_01_04_08.jld"
]

file = files[3]
result = load(file, "result")

#plot(title = "S_pi_0 & S_pi : J = 1; 16 sites lattice")

h_vals = []
m = []
push!(h_vals, result[1]...)
push!(m, result[4]...)

plot!(
h_vals,
m,
dpi = 800,
xlabel = "h",
ylabel = "S_pi",
legend = :topright,
label = "S_pi (L=16)"
)
