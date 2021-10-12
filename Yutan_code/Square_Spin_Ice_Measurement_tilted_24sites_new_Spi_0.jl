using LinearAlgebra
using SparseArrays
using Arpack
using KrylovKit
using Dates
using JLD
using Plots

L1=11;L2=11;J1=1;J2=1;h=0
N = 24

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

function neib(n::Int64;L1::Int64=2, L2::Int64=2)
    coord = coordinate(n,L1=L1, L2=L2)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L2), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L2), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L1)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L1)))
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

function neib_list_gen(;L1::Int64=2, L2::Int64=2)
    neib_list = Set{Int64}[]
    for n in 1:L1*L2
        push!(neib_list, neib(n, L1=L1, L2=L2))
    end
    return neib_list
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

syms = [(3,3),(4,-4)]
# L1 is wide (x direction), L2 is length (y direction)
function translation(coord::Array{Int64, 1}; direction::Symbol, dist::Int64)
    translated_coord = [coord[1], coord[2]]
    if direction == :x
        coord[2] += dist
    end
    if direction == :y
        coord[1] += dist
    end
    return nothing
end

function isvalid(coord;L1::Int64, L2::Int64)
    return (coord[1]>0)&&(coord[2]>0)&&(coord[1]≤L2)&&(coord[2]≤L1)
end

# Represent the symmetry operations by a tuple: (4,2),(-4,4)
function find_translated_coords(coord::Tuple{Int64, Int64},syms::Tuple{Int64, Int64}...;N,L1::Int64, L2::Int64)
    equi_coords = Set()
    push!(equi_coords, coord)
    origin_coord = coord
    coord = [coord[1], coord[2]]
    function loop_over_syms(coord; current_layer = 1, tot_layer = length(syms))
        last_coord = Tuple(coord)
        if current_layer < tot_layer
            loop_over_syms(coord; current_layer = current_layer+1)
            for i in 1:Int(ceil(L1*L2/N))
                sym_op = syms[current_layer]
                translation(coord; direction=:x, dist=i*sym_op[1])
                translation(coord; direction=:y, dist=i*sym_op[2])
                loop_over_syms(coord; current_layer = current_layer+1)
                coord = [last_coord[1], last_coord[2]]

                translation(coord; direction=:x, dist=-i*sym_op[1])
                translation(coord; direction=:y, dist=-i*sym_op[2])
                loop_over_syms(coord; current_layer = current_layer+1)
                coord = [last_coord[1], last_coord[2]]
            end
        else
            if isvalid(coord, L1=L1, L2=L2)
                push!(equi_coords, Tuple(coord))
            end
            coord = [last_coord[1], last_coord[2]]
            for i in 1:Int(ceil(L1*L2/N))
                sym_op = syms[current_layer]
                translation(coord; direction=:x, dist=i*sym_op[1])
                translation(coord; direction=:y, dist=i*sym_op[2])
                if isvalid(coord, L1=L1, L2=L2)
                    push!(equi_coords, Tuple(coord))
                end
                coord = [last_coord[1], last_coord[2]]

                translation(coord; direction=:x, dist=-i*sym_op[1])
                translation(coord; direction=:y, dist=-i*sym_op[2])
                if isvalid(coord, L1=L1, L2=L2)
                    push!(equi_coords, Tuple(coord))
                end
                coord = [last_coord[1], last_coord[2]]
            end
            coord = [origin_coord[1], origin_coord[2]]
        end
    end
    loop_over_syms(coord)
    return equi_coords
end

function find_valid_sites(;L1::Int64, L2::Int64, N::Int64)
    valid_coords = []
    all_coords = []
    for i in 1:L2
        for j in 1:L1
            push!(all_coords, (i, j))
        end
    end
    for coord in all_coords
        i = coord[1]; j = coord[2]
        if (i+j ≥ 10)&&(i+j ≤ 16)&&(i-j ≤ 4)&&(i -j ≥ -4)
            is_valid_site = true
            equi_coords = find_translated_coords(coord,syms...;N=N,L1=L1, L2=L2)
            for equi_coord in equi_coords
                if ((equi_coord[1]+equi_coord[2] ≥ 10)&&(equi_coord[1]+equi_coord[2] ≤16)&&(equi_coord[1]-equi_coord[2] ≤ 4)&&(equi_coord[1] - equi_coord[2] ≥ -4))&&((i,j)>equi_coord)
                    is_valid_site = false
                end
            end
            if is_valid_site
                #println("$(coord) is valid")
                push!(valid_coords, coord)
            end
        else
            continue
        end
    end
    return valid_coords
end

function corre_site_dic(;L1::Int64, L2::Int64, N::Int64, syms, sites)
    site2site = Dict()
    all_coords = []
    for i in 1:L2
        for j in 1:L1
            push!(all_coords, (i, j))
        end
    end
    for coord in all_coords
        i = coord[1]; j = coord[2]

        equi_coords = find_translated_coords(coord,syms...;N=N,L1=L1, L2=L2)
        for equi_coord in equi_coords
            if (equi_coord[1],equi_coord[2]) ∈ sites
                site2site[(i, j)] = equi_coord
            end
        end

    end
    return site2site
end

function tilted_neib_list_gen(;site2site, sites, neib_list, L1::Int64=11, L2::Int64=11)
    new_neib_list = Dict()
    for site in sites
        bit_position = bit_pos((site[1], site[2]);L1=L1,L2=L2)
        neibs = neib_list[bit_position]
        new_neibs = Set()
        for neib in neibs
            if coordinate(neib;L1=L1,L2=L2) ∉ sites
                push!(new_neibs, bit_pos(site2site[coordinate(neib;L1=L1,L2=L2)];L1=L1,L2=L2) )
            else
                push!(new_neibs, neib)
            end
        end
        new_neib_list[bit_pos((site[1], site[2]);L1=L1,L2=L2)] = new_neibs
    end
    return new_neib_list
end

function numbered_tilted_neib_list_gen(tilted_neib_list)
    all_keys = []
    for key in keys(tilted_neib_list)
        push!(all_keys, key)
    end
    sort!(all_keys)

    key2numbering = Dict()
    for i in 1:length(all_keys)
        key2numbering[all_keys[i]] = i
    end

    numbered_tilted_neib_list = []
    for i in 1:length(all_keys)
        neibs = Set{Int64}()
        for value in tilted_neib_list[all_keys[i]]
            push!(neibs, key2numbering[value])
        end
        push!(numbered_tilted_neib_list, neibs)
    end

    return numbered_tilted_neib_list
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

function Hamiltonian1_update(H1;h_new, h_old, nearest_neib_list, second_neib_list)
    h_ratio = h_new/ h_old
    for state in 0:(2^N-1) #loop over all states
        for i in 1:N #loop over all sites in a given state
            flipped_state = state ⊻ (1<<(i-1))
            H1[state+1, flipped_state+1] = H1[state+1, flipped_state+1]*h_ratio
        end
    end
    return H1
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

function Hamiltonian2_update(H2;h_new, h_old, nearest_neib_list, second_neib_list)
    h_ratio = h_new/ h_old
    for state in 0:(2^N-1) #loop over all states
        H2[state+1, state+1] = H2[state+1, state+1] * h_ratio
    end
    return H2
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

function S_pi_H1(state::Array{Float64,1}; N::Int64, neib_list)
    odd_even = [1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,0,1,0,0]

    S_pi = spzeros(2^N,2^N)
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi[basis_state+1,basis_state+1] += (basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*(-1)^(odd_even[i]+odd_even[j])/N
            end
        end
    end
    #now that we have matrix m, calculate the average m_val:
    S_pi_val = conj.(state')*S_pi*state
    return S_pi_val[1] #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_pi_0_H1(state::Array{Float64,1}; N::Int64, neib_list)
    #build up the phase factor for use later
    odd_even = [1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,0,1]
    # start calculating S_pi
    S_pi = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi += conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*(-1)^(odd_even[i]+odd_even[j])/N
            end
        end
    end
    return S_pi #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_0_pi_H1(state::Array{Float64,1}; N::Int64, neib_list)
    #build up the phase factor for use later
    odd_even = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]
    # start calculating S_pi
    S_pi = 0
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            for j in 1:N #loop over all sites again
                 S_pi += conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*(-1)^(odd_even[i]+odd_even[j])/N
            end
        end
    end
    return S_pi #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_real_H1(state::Array{Float64,1}; N::Int64, neib_list)
    # start calculating S_real
    # There are 2 sites per unit cell. So we should calculate them separately.
    S_real1 = zeros(Float64, L2, L1)
    S_real2 = zeros(Float64, L2, L1)
    for basis_state in 0:(2^N-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:N #loop over all sites in a given state
            S_real1[coordinate(i;L1=L1, L2=L2)[1], coordinate(i;L1=L1, L2=L2)[2]] +=  conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[1]-1/2)*(basis_state_binary[i]-1/2)
        end
        for i in 1:N #loop over all sites in a given state
            S_real2[coordinate(i;L1=L1, L2=L2)[1], coordinate(i;L1=L1, L2=L2)[2]] +=  conj(state[basis_state+1])*state[basis_state+1]*(basis_state_binary[3]-1/2)*(basis_state_binary[i]-1/2)
        end
    end
    return (S_real1, S_real2) #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_entangle(state::Array{Float64,1}; N::Int64, neib_list)
    #find the first four spin interacting plaquette
    plaq_sites = [3,4,7,8] #Has to be an array to make it ordered
    env_sites = [1,2,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

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

function driver(;N::Int64=24,L1::Int64, L2::Int64, J1=1, J2, h_vals)
    N = N
    J1 = 1
    #generating neib lists
    syms = [(3,3),(4,-4)]
    nearest_neib_list = nearest_neib_list_gen(L1=L1, L2=L2)
    second_neib_list = second_neib_list_gen(L1=L1, L2=L2)
    all_neib_list = neib_list_gen(L1=L1, L2=L2)
    sites = find_valid_sites(L1=11, L2=11, N=N)
    site2site = corre_site_dic(;L1=11, L2=11, N=N, syms=syms, sites=sites)
    tilted_nearest_neib_list = tilted_neib_list_gen(site2site=site2site, sites=sites, neib_list=nearest_neib_list, L1=11, L2=11)
    tilted_second_neib_list = tilted_neib_list_gen(site2site=site2site, sites=sites, neib_list=second_neib_list, L1=11, L2=11)
    tilted_all_neib_list = tilted_neib_list_gen(site2site=site2site, sites=sites, neib_list=all_neib_list, L1=11, L2=11)
    numbered_tilted_nearest_neib_list = numbered_tilted_neib_list_gen(tilted_nearest_neib_list)
    numbered_tilted_second_neib_list = numbered_tilted_neib_list_gen(tilted_second_neib_list)
    numbered_tilted_all_neib_list = numbered_tilted_neib_list_gen(tilted_all_neib_list)
    #plaquette_list = plaquette_list_gen(;L=L, neib_list=neib_list)
########initilaize containers of observables
    m = []
    S_pi = []
    Fidelity = []
    #Fcl_vals = []
    #Fqm_vals = []
    S_entangle_vals = []
    S_0_pi = []
    S_pi_0 = []
    S_real = []

    eigstate_prev = zeros(Float64, 2^N)
    h_prev = -1
    for i in 1:length(h_vals)
        h = h_vals[i]
        H1 = Hamiltonian1(;N=N, J1=J1, J2=J2, h=h_prev, nearest_neib_list=numbered_tilted_nearest_neib_list, second_neib_list=numbered_tilted_second_neib_list)
        eigstate1 = eigsolve(H1, 1, :SR, eltype(H1), krylovdim = 100)[2][1]
        #calculating Fidelity using eigenstates of H1
        Fid = conj.(eigstate_prev')*eigstate1
        (i == 1) && (Fid[1] = 1) #set Fid to be 1 manually for the first h
        push!(S_pi, S_pi_H1(eigstate1;N=N, neib_list=numbered_tilted_all_neib_list))
        push!(Fidelity, 2*(1-abs(Fid[1]))/(h-h_prev)^2)
        #push!(Fcl_vals, Fcl(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list))
        #push!(Fqm_vals, Fqm(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list))
        #push!(S_entangle_vals, S_entangle(eigstate1;N=N, L1=L1, L2=L2, second_neib_list=second_neib_list))
        eigstate_prev = eigstate1
        push!(S_0_pi, S_0_pi_H1(eigstate1;N=N, neib_list=numbered_tilted_all_neib_list))
        push!(S_pi_0, S_pi_0_H1(eigstate1;N=N, neib_list=numbered_tilted_all_neib_list))
        push!(S_real, S_real_H1(eigstate1; N=N, neib_list=numbered_tilted_all_neib_list))

        #Start Calculating using H2
        #=H2 = Hamiltonian2(;N=N, J1=J1, J2=J2, h=h, nearest_neib_list=nearest_neib_list, second_neib_list=second_neib_list)
        eigstate2 = eigsolve(H2, 1, :SR, eltype(H2), krylovdim = 100)[2][1]
        push!(m, m_H2(eigstate2; N=N))=#
        h_prev = h
    end

    #return (h_vals, m, S_pi, Fidelity, Fcl_vals, Fqm_vals, S_entangle_vals)
    return (h_vals, J2, m, S_pi, Fidelity, S_entangle_vals, S_0_pi, S_pi_0, S_real)
end

L1=11; L2=11; h_vals = range(0.01, 1.5, length = 2); J2_vals = [1.0]
for J2 in J2_vals
    @time result = driver(L1=L1, L2=L2, h_vals=h_vals, J2=J2)

    my_time = Dates.now()

    time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
    content = "Square_Spin_Ice_Measurement_tilted24sites_new"
    save_path = "E:/UC Davis/Research/Square Spin Ice/Square-Spin-Ice/Yutan_code/Results/"
    #"/nfs/home/zyt329/Research/Square_spin_ice/result/"
    save_name = save_path*content*"_J2=$(J2)_hmin=$(h_vals[1])_hmax=$(h_vals[end])_"*time_finished*".jld"

    save(save_name, "result", result)
    println("J2 = $(J2) finished")
end
