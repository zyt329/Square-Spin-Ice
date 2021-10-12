using JLD
using Plots

save_path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/plots/"

files = ["E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_J2=1.0_L1=4_L2=6_hmin=0.01_hmax=1.5_Date_Mon_04_Oct_2021_13_24_14.jld"
]

file = files[1]
result = load(file, "result")

plot(title = "S_pi_0 & S_pi : J = 1.0; 24 sites lattice")

h_vals = []
m = []
push!(h_vals, result[1]...)
push!(m, result[9]...)

plot!(
h_vals,
m,
dpi = 800,
xlabel = "h",
ylabel = "S_pi_0",
legend = :topright,
label = "S_pi_0 (L=24)"
)
