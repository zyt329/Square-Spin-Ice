using JLD
using Dates
using LinearAlgebra
using SparseArrays
using Arpack

L=4

function coordinate(n;L::Int64=2)
    num_sites = L^2
    i::Int64 = Int(ceil(n/L))
    j::Int64 = mod1(n,L)  #site i is at i-th row, j-th column
    return (i,j)
end

function bit_pos(coordinate::Tuple{Int64,Int64};L::Int64=2)
    n = (coordinate[1]-1)*L + coordinate[2]
    return n
end

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

neib_list = neib_list_gen(L = L)

function translation(state::Int64; direction::Symbol, dist::Int64, L::Int64)
    max_len = L^2
    state_binary = digits!(zeros(Int64, 64), state, base = 2)
    translated_state = 0
    if direction == :x
        for pos in 1:max_len
            if state_binary[pos] == 1
                #translated_pos = mod1(pos + dist, L) + Int(floor((pos-0.5)/L))*L
                #translated_state = translated_state ⊻ (1<<(translated_pos-1))
                pos_coord = coordinate(pos,L=L)
                trans_pos_coord_i = pos_coord[1]
                trans_pos_coord_j = mod1(pos_coord[2] + dist, L)
                trans_pos = bit_pos((trans_pos_coord_i, trans_pos_coord_j), L=L)
                translated_state = translated_state ⊻ (1<<(trans_pos-1))
            end
        end
        return translated_state
    elseif direction == :y
        for pos in 1:max_len
            if state_binary[pos] == 1
                #translated_pos = Int((mod1(ceil((pos-0.5)/L) + dist, L)-1)*L + mod1(pos, L))
                #translated_state = translated_state ⊻ (1<<(translated_pos-1))
                pos_coord = coordinate(pos,L=L)
                trans_pos_coord_i = mod1(pos_coord[1] + dist, L)
                trans_pos_coord_j = pos_coord[2]
                trans_pos = bit_pos((trans_pos_coord_i, trans_pos_coord_j), L=L)
                translated_state = translated_state ⊻ (1<<(trans_pos-1))
            end
        end
        return translated_state
    else
        error("Direction not defined.")
    end
end

function checkstate(s::Int64; k::Real, direction::Symbol, L::Int64)
    R = -1; t = s;
    for i in 1:L
        t = translation(t, direction=direction, dist=1; L=L)
        if t < s
            #println("t<s")
            return (false, R)
        elseif t == s
            R = i
            if mod(k, L/i) != 0
                #println("$(L/i)")
                return (false, R)
            else
                #println("should be true")
                return (true, R)
            end
        end
    end
end

#Check whether a state is ref state 2D, if it's a ref state, return number of different states Da
function checkstate(s::Int64; k::Tuple, L::Int64, ϵ::Float64 = 10^(-10))
    k = 2pi/L .* k
    states = Set{Int64}()
    push!(states, s)
    is_ref = true
    Fa = 0
    t = s
    for i in 1:L
        t = translation(t, direction=:x, dist=1; L=L)
        t = translation(t, direction=:y, dist=1; L=L)
        for j in 1:L
            t = translation(t, direction=:x, dist=-1; L=L)
            t = translation(t, direction=:y, dist=1; L=L)
            push!(states, t)
            # Check if the state is the smallest integer and return periodicity Rx, Ry
            if t < s
                is_ref = false
            end
            if t == s
                Fa += exp(-im*((k[1] + k[2])*i + (-k[1] + k[2])*j))
            end
            #finished check states for one certain translation
        end
    end
    if abs(Fa) < ϵ
       is_ref = false
    end
    return (is_ref, length(states), real(Fa))
end

# generate a list of ref states for a 2D system, with k being a 2D Tuple (kx, ky)
# return (reference state, number of different states for this reference state, sum of phase for the k state)
function ref_state_gen(;k::Tuple, L::Int64)
    ref_states = Int64[]#Set{Int64}()
    Da = Int64[]
    Fa = Number[]
    ref_state_nums = Dict{Int64, Int64}()
    ref_state_tots = 0
    for state in 0:(2^(L^2)-1)
        chk_state = checkstate(state, k=k, L=L)
        if chk_state[1]
            ref_state_tots += 1
            ref_state_nums[state] = ref_state_tots
            push!(ref_states, state)
            push!(Da, chk_state[2])
            push!(Fa, chk_state[3])
        end
    end
    return Dict(:ref_states => ref_states, :Da => Da, :Fa => Fa, :ref_state_tots => ref_state_tots, :ref_state_nums => ref_state_nums)
end

# return the number of translation in x and y direction in a Tuple (i, j)
function find_ref(s::Int64; L::Int64)
    t = s
    ref = t
    trans = (0,0)
    for i in 1:L
        t = translation(t, direction=:x, dist=1; L=L)
        t = translation(t, direction=:y, dist=1; L=L)
        for j in 1:L
            t = translation(t, direction=:x, dist=-1; L=L)
            t = translation(t, direction=:y, dist=1; L=L)
            if t < ref
                ref = t
                trans = (i-j, i+j)
            end
        end
    end
    return Dict(:trans => trans, :ref => ref)
end

function Hamiltonian(;k::Tuple{Int64, Int64},L::Int64,J=1, h=1, neib_list, state_list, state_nums, state_tot, Das, Fas)
    H = zeros(Complex{Float64}, state_tot,state_tot)
    for state in state_list #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        Na = Das[state_nums[state]]*abs(Fas[state_nums[state]])^2 # for calculation of Hamiltonian later
        for i in 1:L^2 #loop over all sites i in a given state
            #diagonal terms interact with h
            if state_binary[i] == 1
                H[state_nums[state],state_nums[state]] += h
            else
                H[state_nums[state],state_nums[state]] -= h
            end
            # off diagonal terms come from flipping bonds, j is neighbor of i
            for j in neib_list[i]
                flipped_state = state ⊻ (1<<(i-1))
                flipped_state = flipped_state ⊻ (1<<(j-1))
                ref_flipped_state = find_ref(flipped_state, L = L)
                if haskey(state_nums, ref_flipped_state[:ref]) # check if the flipped state is within this k sector
                    # calculated the contribution to the Hamiltonian
                    Nb = Das[state_nums[ref_flipped_state[:ref]]]*abs(Fas[state_nums[ref_flipped_state[:ref]]])^2
                    H[state_nums[state],state_nums[ref_flipped_state[:ref]]] += 1/2 * J/4 * exp(2pi*im*(k[1]*ref_flipped_state[:trans][1] + k[2]*ref_flipped_state[:trans][2]) / L) * sqrt(Nb / Na)
                end
            end
        end
    end
    return H
end

function hermitian_test()
    eigs = []
    for kx in 1:L
        for ky in 1:L
            ref_state = ref_state_gen(k=(kx,ky),L=L)
            H = Hamiltonian(;k=(kx,ky),L=L,J=1, h=0, neib_list=neib_list, state_list=ref_state[:ref_states], state_nums=ref_state[:ref_state_nums], state_tot=ref_state[:ref_state_tots], Das=ref_state[:Da], Fas=ref_state[:Fa])
            if maximum(abs.(Hermitian(H, :U) - Hermitian(H, :L))) > 10^(-10)
                return false
            end
            #push!(eigs, eigen(Hermitian(H,:U)).values...)
        end
    end
    return true
end

allowed_k = [(0,0),(0,1),(1,0),(-1,0),(0,-1),(2,0),(1,1),(-1,1)]

function Hamiltonian_solver(allowed_k)
    eigs = []
    for k in allowed_k
        println("Dealing with k = ($(k[1]),$(k[2]))")
        ref_state = ref_state_gen(k=k,L=L)
        H = Hamiltonian(;k=k,L=L,J=1, h=0, neib_list=neib_list, state_list=ref_state[:ref_states], state_nums=ref_state[:ref_state_nums], state_tot=ref_state[:ref_state_tots], Das=ref_state[:Da], Fas=ref_state[:Fa])
        if maximum(abs.(Hermitian(H, :U) - Hermitian(H, :L))) > 10^(-10)
            return ((k[1],k[2]),false, maximum(abs.(Hermitian(H, :U) - Hermitian(H, :L))))
        end
        push!(eigs, eigen(Hermitian(H,:U)))
    end
    return eigs
end

eigens = Hamiltonian_solver(allowed_k)

my_time = Dates.now()

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
content = "Momentum_Block_Diag_frustrated_lattice"
save_path = "/nfs/home/zyt329/Research/Square_spin_ice/result/"
#"D:/UC Davis/Research/Square Spin Ice/Square-Spin-Ice/Yutan_code/Results/"
save_name = save_path*content*"_L=$(L)__J=1__h=0_"*time_finished*".jld"

save(save_name, "eigens", eigens)
println("finished")
