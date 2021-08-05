using LinearAlgebra
using SparseArrays
using Arpack
using KrylovKit
using Dates
using JLD

function coordinate(n;L::Int64=2)
    num_sites = L^2
    i::Int64 = Int(ceil(n/L))
    j::Int64 = mod1(n,L)  #site i is at i-th row, j-th column
    return (i,j)
end

function bit_pos(coordinate::Tuple{Int64,Int64};L::Int64=2)
    n = (coordinate[1]-1)*L + coordinate[2]
    return n
end# Now Let's Construct Hamiltonian for our problem, by changing neighbor lists

function neib(n::Int64;L::Int64=2)
    coord = coordinate(n,L=L)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L)))
    if iseven(coord[1]+coord[2])
        push!(neibs, (mod1(coord[1]+1,L), mod1(coord[2]-1,L)))
        push!(neibs, (mod1(coord[1]-1,L), mod1(coord[2]+1,L)))
    else
        push!(neibs, (mod1(coord[1]+1,L), mod1(coord[2]+1,L)))
        push!(neibs, (mod1(coord[1]-1,L), mod1(coord[2]-1,L)))
    end
    #=convert coordinations to positions in bits=#
    neibs_bit_pos = Set{Int64}()
    for neib in neibs
        push!(neibs_bit_pos, bit_pos(neib,L=L))
    end
    return neibs_bit_pos
end

function neib_list_gen(;L::Int64=L)
    neib_list = Set{Int64}[]
    for n in 1:L^2
        push!(neib_list, neib(n, L=L))
    end
    return neib_list
end

function Hamiltonian1(;L::Int64=2, J=1, h=1, neib_list)
    H = spzeros(2^(L^2),2^(L^2))
    for state in 0:(2^(L^2)-1) #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for i in 1:L^2 #loop over all sites in a given state
            flipped_state = state ⊻ (1<<(i-1))
            H[state+1,flipped_state+1] += (1/2)*h
            for j in neib_list[i] #loop over(compare) all neighbors of a given site
                H[state+1,state+1] += (state_binary[i]-1/2)*(state_binary[j]-1/2)*J/2
            end
        end
    end
    return H
end

function Hamiltonian2(;L::Int64=2, J=1, h=1, neib_list)
    H = spzeros(2^(L^2),2^(L^2))
    for state in 0:(2^(L^2)-1) #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for i in 1:L^2 #loop over all sites in a given state
            if state_binary[i] == 1
                H[state+1,state+1] -= h/2
            else
                H[state+1,state+1] += h/2
            end
            for j in neib_list[i] #loop over(compare) all neighbors of a given site
                flipped_state = state ⊻ (1<<(i-1))
                flipped_state = flipped_state ⊻ (1<<(j-1))
                H[state+1,flipped_state+1] += (1/4)*(J/2)
            end
        end
    end
    return H
end

function driver(;L=L, h_vals)
    J = 1
    neib_list = neib_list_gen(L=L)
    H1_eigenstates = []
    H2_eigenstates = []
    for i in 1:length(h_vals)
        h = h_vals[i]
        H1 = Hamiltonian1(;L=L, J=J, h=h, neib_list=neib_list)
        H2 = Hamiltonian2(;L=L, J=J, h=h, neib_list=neib_list)
        push!(H1_eigenstates, eigsolve(H1, 1, :SR, eltype(H1)))
        push!(H2_eigenstates, eigsolve(H2, 1, :SR, eltype(H2)))
    end

    return (H1_eigenstates, H2_eigenstates, h_vals)
end

L=4; h_vals = range(0.001, 1, length = 200)
@time result = driver(L=L, h_vals=h_vals)

my_time = Dates.now()

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
content = "Ground_state_generator"
save_path = "/nfs/home/zyt329/Research/Square_spin_ice/result/"
#"E:/UC Davis/Research/Square Spin Ice/Square-Spin-Ice/Yutan_code/Results/"
save_name = save_path*content*"_L=$(L)_hmax=$(h_vals[end])_"*time_finished*".jld"

save(save_name, "result", result)
println("finished")
