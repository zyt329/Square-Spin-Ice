using LinearAlgebra
using SparseArrays
using Arpack
using JLD
using Dates
using KrylovKit

simulation = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_J2=1_L1=4_L2=6_hmin=0.01_hmax=100.0_Date_Sun_08_Aug_2021_18_40_11.jld"

result = load(simulation, "result")
#=h_vals = result[1]
m = result[2]
S_pi = result[3]
Fidelity = result[4]
Fcl_vals = result[5]
Fqm_vals = result[6]=#
S_entangle_vals = result[6]

println(S_entangle_vals[1])
#=H_even = load(simulation)["sim"][1]
H_odd = load(simulation)["sim"][2]

println(eigsolve(H_odd, 3, :SR, krylovdim=200)[1])=#
