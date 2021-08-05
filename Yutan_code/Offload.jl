using LinearAlgebra
using SparseArrays
using Arpack
using JLD
using Dates
using KrylovKit

simulation = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Square_Spin_Ice_Measurement_L=4_hmin=99.99999999999999_hmax=1.0e8_Date_Tue_27_Jul_2021_16_34_33.jld"

result = load(simulation, "result")
h_vals = result[1]
m = result[2]
S_pi = result[3]
Fidelity = result[4]
Fcl_vals = result[5]
Fqm_vals = result[6]
S_entangle_vals = result[7]

println(S_pi[end])
#=H_even = load(simulation)["sim"][1]
H_odd = load(simulation)["sim"][2]

println(eigsolve(H_odd, 3, :SR, krylovdim=200)[1])=#
