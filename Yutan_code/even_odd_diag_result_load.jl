using JLD
using Plots

path = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/"
name = "Hamiltonian_even_odd_L=4__J=1.0__h=1.0_Date_Thu_17_Jun_2021_22_39_55.jld"
result = load(path*name, "eigs")
