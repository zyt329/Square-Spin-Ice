using LinearAlgebra
using SparseArrays
using Arpack
using KrylovKit
using Dates
using JLD

function coordinate(n;L1::Int64=2,L2::Int64=2)
    @assert((n ≤ L1 * L2) && (1 ≤ n),"The numbering (bit position) of a site shouldn't exceed the total number of sites $(L1 * L2), and should be bigger than 0.")
    i::Int64 = Int(ceil(n/L1))
    j::Int64 = mod1(n,L1)  #site i is at i-th row, j-th column
    return (i,j)
end

function bit_pos(coordinate::Tuple{Int64,Int64};L1::Int64=2,L2::Int64=2)
    @assert((coordinate[1] ≤ L2) && (coordinate[2] ≤ L1),"The cooridnate should be within the range of the lattice size $L1 by $L2")
    n = (coordinate[1]-1)*L1 + coordinate[2]
    return n
end

function nearest_neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L2), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L2), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L1)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L1)))
    #=convert coordinations to positions in bits=#
    neibs_bit_pos = Set{Int64}()
    for neib in neibs
        push!(neibs_bit_pos, bit_pos(neib, L1=L1, L2=L2))
    end
    return neibs_bit_pos
end

function second_neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    if iseven(coord[1]+coord[2])
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]-1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]+1,L1)))
    else
        push!(neibs, (mod1(coord[1]+1,L2), mod1(coord[2]+1,L1)))
        push!(neibs, (mod1(coord[1]-1,L2), mod1(coord[2]-1,L1)))
    end
    #=convert coordinations to positions in bits=#
    neibs_bit_pos = Set{Int64}()
    for neib in neibs
        push!(neibs_bit_pos, bit_pos(neib, L1=L1, L2=L2))
    end
    return neibs_bit_pos
end

function nearest_neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, nearest_neib(n, L1=L1, L2=L2))
    end
    return neib_list
end

function second_neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, second_neib(n, L1=L1, L2=L2))
    end
    return neib_list
end

function chk_parity(state::Int64)
    state_binary = digits!(zeros(Int64, 64), state, base = 2)
    if iseven(sum(state_binary))
        return :even
    else
        return :odd
    end
end

function parity_states_list(;N::Int64, L1::Int64=L1, L2::Int64=L2)
    even_state = Int64[]
    odd_state = Int64[]
    even_state_num = Dict{Int64, Int64}()
    odd_state_num = Dict{Int64, Int64}()
    even_state_tot = 0
    odd_state_tot = 0
    for state in 0:(2^N-1)
        if chk_parity(state) == :even
            even_state_tot += 1
            push!(even_state, state)
            even_state_num[state] = even_state_tot
        else
            odd_state_tot += 1
            push!(odd_state, state)
            odd_state_num[state] = odd_state_tot
        end
    end
    return Dict{Symbol, Any}(:even_state => even_state, :odd_state => odd_state, :even_state_num => even_state_num, :odd_state_num => odd_state_num, :even_state_tot => even_state_tot, :odd_state_tot => odd_state_tot)
end

function update_val(row_inds, col_inds, vals;row_ind, col_ind, val)
    push!(row_inds, row_ind)
    push!(col_inds, col_ind)
    push!(vals, val)
end

function Hamiltonian1(;N::Int64=2, J1=1, J2=1, h=1, nearest_neib_list, second_neib_list)
    row_inds = Int64[]
    col_inds = Int64[]
    vals = Float64[]
    for state in 0:(2^N-1) #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for i in 1:N #loop over all sites in a given state
            flipped_state = state ⊻ (1<<(i-1))
            update_val(row_inds, col_inds, vals, row_ind=state+1, col_ind = flipped_state+1, val = (1/2)*h)
            for j in nearest_neib_list[i] #loop over(compare) all neighbors of a given site
                update_val(row_inds, col_inds, vals, row_ind = state+1, col_ind = state+1, val = (state_binary[i]-1/2)*(state_binary[j]-1/2)*J1/2)
            end
            for j in second_neib_list[i] #loop over(compare) all neighbors of a given site
                update_val(row_inds, col_inds, vals, row_ind = state+1, col_ind = state+1, val = (state_binary[i]-1/2)*(state_binary[j]-1/2)*J2/2)
            end
        end
    end
    return sparse(row_inds, col_inds, vals, 2^N, 2^N, +)
end

function Hamiltonian2(;N::Int64=2, J1=1, J2=1, h=1, nearest_neib_list, second_neib_list)
    row_inds = Int64[]
    col_inds = Int64[]
    vals = Float64[]
    for state in 0:(2^N-1) #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for i in 1:N #loop over all sites in a given state
            if state_binary[i] == 1
                update_val(row_inds, col_inds, vals, row_ind=state+1, col_ind = state+1, val = - h/2)
            else
                update_val(row_inds, col_inds, vals, row_ind=state+1, col_ind = state+1, val =  h/2)
            end
            for j in nearest_neib_list[i] #loop over(compare) all neighbors of a given site
                flipped_state = state ⊻ (1<<(i-1))
                flipped_state = flipped_state ⊻ (1<<(j-1))
                update_val(row_inds, col_inds, vals, row_ind=state+1, col_ind = flipped_state+1, val =  (1/4)*(J1/2))
            end
            for j in second_neib_list[i] #loop over(compare) all neighbors of a given site
                flipped_state = state ⊻ (1<<(i-1))
                flipped_state = flipped_state ⊻ (1<<(j-1))
                update_val(row_inds, col_inds, vals, row_ind=state+1, col_ind = flipped_state+1, val =  (1/4)*(J2/2))
            end
        end
    end
    return sparse(row_inds, col_inds, vals, 2^N, 2^N, +)
end

function m_H2(state::Array{Float64,1}; N::Int64)
    m = spzeros(2^N,2^N)
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        # calculating total spin along z direction, considering it's spin 1/2
        m[basis_state+1, basis_state+1] += (sum(basis_state_binary)-1/2*N)/N
    end
    #now that we have matrix m, calculate the average m_val:
    m_val = conj.(state')*m*state
    return m_val[1] #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function m_H2(state::Array{Float64,1}; N::Int64)
    m = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        # calculating total spin along z direction, considering it's spin 1/2
        m += conj(state[basis_state+1])*state[basis_state+1]*(sum(basis_state_binary)-1/2*N)/N
    end
    return m
end

function S_pi_H1(state::Array{Float64,1}; N::Int64, L1::Int64, L2::Int64)
    S_pi = spzeros(2^N,2^N)
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi[basis_state+1,basis_state+1] += (basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*(-1)^(sum(coordinate(i;L1=L1, L2=L2))+sum(coordinate(j;L1=L1, L2=L2)))/N
            end
        end
    end
    #now that we have matrix m, calculate the average m_val:
    S_pi_val = conj.(state')*S_pi*state
    return S_pi_val[1] #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_pi_H1(state::Array{Float64,1}; N::Int64, L1::Int64, L2::Int64)
    #build up the phase factor for use later
    phase_factor = zeros(Int64, N, N)
    for i in 1:N #loop over all sites in a given state
        for j in 1:N #loop over all sites again
            phase_factor[i, j] = (-1)^(sum(coordinate(i;L1=L1, L2=L2))+sum(coordinate(j;L1=L1, L2=L2)))
        end
    end
    # start calculating S_pi
    S_pi = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi += conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*phase_factor[i, j]/N
            end
        end
    end
    return S_pi #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_pi_0_H1(state::Array{Float64,1}; N::Int64, L1::Int64, L2::Int64)
    #build up the phase factor for use later
    phase_factor = zeros(Int64, N, N)
    for i in 1:N #loop over all sites in a given state
        for j in 1:N #loop over all sites again
            phase_factor[i, j] = (-1)^(coordinate(i;L1=L1, L2=L2)[2]-coordinate(j;L1=L1, L2=L2)[2])
        end
    end
    # start calculating S_pi
    S_pi = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi += conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*phase_factor[i, j]/N
            end
        end
    end
    return S_pi #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_0_pi_H1(state::Array{Float64,1}; N::Int64, L1::Int64, L2::Int64)
    #build up the phase factor for use later
    phase_factor = zeros(Int64, N, N)
    for i in 1:N #loop over all sites in a given state
        for j in 1:N #loop over all sites again
            phase_factor[i, j] = (-1)^(coordinate(i;L1=L1, L2=L2)[1]-coordinate(j;L1=L1, L2=L2)[1])
        end
    end
    # start calculating S_pi
    S_pi = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi += conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*phase_factor[i, j]/N
            end
        end
    end
    return S_pi #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_entangle(state::Array{Float64,1}; N::Int64, L1::Int64, L2::Int64, second_neib_list)
    #find the first four spin interacting plaquette
    plaq_sites = Int64[] #Has to be an array to make it ordered
    for site in 1:N #loop over all sites
        site_coord = coordinate(site,L1=L1, L2=L2)
        if bit_pos((mod1(site_coord[1]-1, L2),mod1(site_coord[2]-1, L1)), L1=L1, L2=L2) ∈ second_neib_list[site]
            push!(plaq_sites, site)
            push!(plaq_sites, bit_pos((mod1(site_coord[1]+1, L2),site_coord[2]), L1=L1, L2=L2))
            push!(plaq_sites, bit_pos((mod1(site_coord[1]+1, L2),mod1(site_coord[2]+1, L1)), L1=L1, L2=L2))
            push!(plaq_sites, bit_pos((site_coord[1],mod1(site_coord[2]+1, L1)), L1=L1, L2=L2))
            break
        end
    end
    #println(plaq_sites)
    env_sites = Int64[]
    for site in 1:N
        if site ∉ plaq_sites
            push!(env_sites, site)
        end
    end
    #println(env_sites)
    C = zeros(Number, 16, Int(2^N/16))
    for basis_state in 0:(2^N - 1)
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        a = 1; b = 1 # start at 1 to avoid 0 as numeration number
        for i in 1:4
            a += basis_state_binary[plaq_sites[i]] * 2^(i-1)
        end
        for i in 1:(N-4)
            b += basis_state_binary[env_sites[i]] * 2^(i-1)
        end
        C[a, b] = state[basis_state + 1]
    end
    Sing_vals = svd(C).S
    return  -sum(Sing_vals.^2 .* log.(Sing_vals.^2))
end

function driver(;L1::Int64, L2::Int64, h_vals, J2)
    N = L1*L2
    J1 = 1
    nearest_neib_list = nearest_neib_list_gen(L1=L1, L2=L2)
    second_neib_list = second_neib_list_gen(L1=L1, L2=L2)
    #plaquette_list = plaquette_list_gen(;L=L, neib_list=neib_list)
    m = []
    S_pi = []
    Fidelity = []
    #Fcl_vals = []
    #Fqm_vals = []
    S_entangle_vals = []
    S_0_pi = []
    S_pi_0 = []

    eigstate_prev = zeros(Float64, 2^N)
    h_prev = -1
    for i in 1:length(h_vals)
        h = h_vals[i]
        H1 = Hamiltonian1(;N=N, J1=J1, J2=J2, h=h, nearest_neib_list=nearest_neib_list, second_neib_list=second_neib_list)
        eigstate1 = eigsolve(H1, 1, :SR, eltype(H1), krylovdim = 100)[2][1]
        #calculating Fidelity using eigenstates of H1
        Fid = conj.(eigstate_prev')*eigstate1
        (i == 1) && (Fid[1] = 1) #set Fid to be 1 manually for the first h
        push!(S_pi, S_pi_H1(eigstate1;N=N, L1=L1, L2=L2))
        push!(Fidelity, 2*(1-abs(Fid[1]))/(h-h_prev)^2)
        #push!(Fcl_vals, Fcl(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list))
        #push!(Fqm_vals, Fqm(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list))
        #push!(S_entangle_vals, S_entangle(eigstate1;N=N, L1=L1, L2=L2, second_neib_list=second_neib_list))
        eigstate_prev = eigstate1
        push!(S_0_pi, S_0_pi_H1(eigstate1;N=N, L1=L1, L2=L2))
        push!(S_pi_0, S_pi_0_H1(eigstate1;N=N, L1=L1, L2=L2))

        #Start Calculating using H2
        #=H2 = Hamiltonian2(;N=N, J1=J1, J2=J2, h=h, nearest_neib_list=nearest_neib_list, second_neib_list=second_neib_list)
        eigstate2 = eigsolve(H2, 1, :SR, eltype(H2), krylovdim = 100)[2][1]
        push!(m, m_H2(eigstate2; N=N))=#
        h_prev = h
    end
    #return (h_vals, m, S_pi, Fidelity, Fcl_vals, Fqm_vals, S_entangle_vals)
    return (h_vals, J2, m, S_pi, Fidelity, S_entangle_vals, S_0_pi, S_pi_0)
end

L1=4; L2=4; h_vals = range(0.01, 10, length = 40); J2_vals = [1.0]
for J2 in J2_vals
    @time result = driver(L1=L1, L2=L2, h_vals=h_vals, J2=J2)

    my_time = Dates.now()

    time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
    content = "Square_Spin_Ice_Measurement"
    save_path = "E:/UC Davis/Research/Square Spin Ice/Square-Spin-Ice/Yutan_code/Results/"
    #"/nfs/home/zyt329/Research/Square_spin_ice/result/"
    save_name = save_path*content*"_J2=$(J2)_L1=$(L1)_L2=$(L2)_hmin=$(h_vals[1])_hmax=$(h_vals[end])_"*time_finished*".jld"

    #save(save_name, "result", result)
    println("J2 = $(J2) finished")
end
