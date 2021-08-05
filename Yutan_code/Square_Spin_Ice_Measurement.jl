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

function m_H2(state::Array{Float64,1}; L::Int64)
    m = spzeros(2^(L^2),2^(L^2))
    for basis_state in 0:(2^(L^2)-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        # calculating total spin along z direction, considering it's spin 1/2
        m[basis_state+1, basis_state+1] += (sum(basis_state_binary)-1/2*L^2)/L^2
    end
    #now that we have matrix m, calculate the average m_val:
    m_val = conj.(state')*m*state
    return m_val[1] #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function S_pi_H1(state::Array{Float64,1}; L::Int64, neib_list)
    S_pi = spzeros(2^(L^2),2^(L^2))
    for basis_state in 0:(2^(L^2)-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        for i in 1:L^2 #loop over all sites in a given state
            for j in 1:L^2 #loop over all sites again
                 S_pi[basis_state+1,basis_state+1] += (basis_state_binary[i]-1/2)*(basis_state_binary[j]-1/2)*(-1)^(sum(coordinate(i;L=L))+sum(coordinate(j;L=L)))/(L^2)
            end
        end
    end
    #now that we have matrix m, calculate the average m_val:
    S_pi_val = conj.(state')*S_pi*state
    return S_pi_val[1] #taking the 1st value of m_val because it's recognized as a length 1 Array
end

function plaquette_list_gen(;L, neib_list)
    plaquette_list = []
    for site in 1:L^2 #loop over all sites
        site_coord = coordinate(site,L=L)
        if bit_pos((mod1(site_coord[1]-1, L),mod1(site_coord[2]-1, L)),L=L) ∈ neib_list[site]
            plaq_site = Int64[] #Has to be an array to make it ordered
            push!(plaq_site, site)
            push!(plaq_site, bit_pos((site_coord[1],mod1(site_coord[2]-1, L)),L=L))
            push!(plaq_site, bit_pos((mod1(site_coord[1]+1, L),mod1(site_coord[2]-1, L)),L=L))
            push!(plaq_site, bit_pos((mod1(site_coord[1]+1, L),site_coord[2]),L=L))
            #push the found plaquette to the plaquette list
            push!(plaquette_list, plaq_site)
        end
    end
    return plaquette_list
end

function plaquette_phase_factors_gen(;L::Int64, plaquette_list)
    @assert(iseven(L), "L must be even for our lattice")
    phase_factors = Int64[]
    for p in 1:Int(L/2)
        for j in 1:Int(L/2)
            push!(phase_factors, 0)
        end
        for k in 1:Int(L/2)
            push!(phase_factors, 1)
        end
    end
    return phase_factors
end

function is_flippable(state_binary::Array, plaquette_sites::Array{Int64, 1})
    spin_plaquette = Int64[]
    for site in plaquette_sites
        push!(spin_plaquette, state_binary[site])
    end
    if (spin_plaquette - [1;0;1;0] == [0;0;0;0]) || (spin_plaquette - [0;1;0;1] == [0;0;0;0])
        return true
    else
        return false
    end
end

function Fcl(state::Array{Float64,1}; L::Int64, neib_list, plaquette_list, plaquette_phase_factors)
    Fcl = spzeros(2^(L^2),2^(L^2))
    for basis_state in 0:(2^(L^2)-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        plaquette_list_len = length(plaquette_list)
        for i in 1:plaquette_list_len
            for j in 1:plaquette_list_len
                if is_flippable(basis_state_binary, plaquette_list[i]) && is_flippable(basis_state_binary, plaquette_list[j])
                    Fcl[basis_state+1, basis_state+1] += (-1)^(plaquette_phase_factors[i]+plaquette_phase_factors[j]) / plaquette_list_len
                end
            end
        end
    end
    # calculating Fcl value for the state
    Fcl_val = conj.(state')*Fcl*state
    return Fcl_val[1]
end

function Fqm(state::Array{Float64,1}; L::Int64, neib_list, plaquette_list, plaquette_phase_factors)
    Fqm = spzeros(2^(L^2),2^(L^2))
    for basis_state in 0:(2^(L^2)-1) #loop over all states
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        plaquette_list_len = length(plaquette_list)
        for i in 1:plaquette_list_len
            for j in 1:plaquette_list_len
                if is_flippable(basis_state_binary, plaquette_list[i]) && is_flippable(basis_state_binary, plaquette_list[j])
                    flipped_state = basis_state
                    for site in vcat(plaquette_list[i], plaquette_list[j])
                        flipped_state = flipped_state ⊻ (1<<(site-1))
                    end
                    Fqm[basis_state+1, flipped_state+1] += (-1)^(plaquette_phase_factors[i]+plaquette_phase_factors[j]) / plaquette_list_len
                end
            end
        end
    end
    # calculating Fcl value for the state
    Fqm_val = conj.(state')*Fqm*state
    return Fqm_val[1]
end

function S_entangle(state::Array{Float64,1}; L::Int64, neib_list)
    #find the first four spin interacting plaquette
    plaq_sites = Int64[] #Has to be an array to make it ordered
    for site in 1:L^2 #loop over all sites
        site_coord = coordinate(site,L=L)
        if bit_pos((mod1(site_coord[1]-1, L),mod1(site_coord[2]-1, L)),L=L) ∈ neib_list[site]
            push!(plaq_sites, site)
            push!(plaq_sites, bit_pos((mod1(site_coord[1]+1, L),site_coord[2]),L=L))
            push!(plaq_sites, bit_pos((mod1(site_coord[1]+1, L),mod1(site_coord[2]+1, L)),L=L))
            push!(plaq_sites, bit_pos((site_coord[1],mod1(site_coord[2]+1, L)),L=L))
            break
        end
    end
    #println(plaq_sites)
    env_sites = Int64[]
    for site in 1:L^2
        if site ∉ plaq_sites
            push!(env_sites, site)
        end
    end
    #println(env_sites)
    C = zeros(Number, 16, Int(2^(L^2)/16))
    for basis_state in 0:(2^(L^2) - 1)
        basis_state_binary = digits!(zeros(Int64, 64), basis_state, base = 2)
        a = 1; b = 1 # start at 1 to avoid 0 as numeration number
        for i in 1:4
            a += basis_state_binary[plaq_sites[i]] * 2^(i-1)
        end
        for i in 1:(L^2-4)
            b += basis_state_binary[env_sites[i]] * 2^(i-1)
        end
        C[a, b] = state[basis_state + 1]
    end
    Sing_vals = svd(C).S
    return  -sum(Sing_vals.^2 .* log.(Sing_vals.^2))
end

function driver(;L=L, h_vals)
    J = 1
    neib_list = neib_list_gen(L=L)
    plaquette_list = plaquette_list_gen(;L=L, neib_list=neib_list)
    plaquette_phase_factors = plaquette_phase_factors_gen(;L=L, plaquette_list = plaquette_list)
    m = []
    S_pi = []
    Fidelity = []
    Fcl_vals = []
    Fqm_vals = []
    S_entangle_vals = []
    eigstate_prev = zeros(Float64, 2^(L^2))
    h_prev = -1
    for i in 1:length(h_vals)
        h = h_vals[i]
        H1 = Hamiltonian1(;L=L, J=J, h=h, neib_list=neib_list)
        H2 = Hamiltonian2(;L=L, J=J, h=h, neib_list=neib_list)
        eigstate1 = eigsolve(H1, 1, :SR, eltype(H1))[2][1]
        eigstate2 = eigsolve(H2, 1, :SR, eltype(H2))[2][1]
        #calculating Fidelity using eigenstates of H1
        Fid = conj.(eigstate_prev')*eigstate1
        (i == 1) && (Fid[1] = 1) #set Fid to be 1 manually for the first h
        push!(m, m_H2(eigstate2; L=L))
        push!(S_pi, S_pi_H1(eigstate1; L=L, neib_list=neib_list))
        push!(Fidelity, (1-abs(Fid[1]))/(h-h_prev)^2)
        push!(Fcl_vals, Fcl(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list,plaquette_phase_factors=plaquette_phase_factors))
        push!(Fqm_vals, Fqm(eigstate1; L=L, neib_list=neib_list, plaquette_list=plaquette_list,plaquette_phase_factors=plaquette_phase_factors))
        push!(S_entangle_vals, S_entangle(eigstate1; L=L, neib_list=neib_list))
        eigstate_prev = eigstate1
        h_prev = h
    end
    return (h_vals, m, S_pi, Fidelity, Fcl_vals, Fqm_vals, S_entangle_vals)
end

L=4; h_vals = range(1, 100, length = 50)
@time result = driver(L=L, h_vals=h_vals)

my_time = Dates.now()

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
content = "Square_Spin_Ice_Measurement"
save_path = "/nfs/home/zyt329/Research/Square_spin_ice/result/"
#"E:/UC Davis/Research/Square Spin Ice/Square-Spin-Ice/Yutan_code/Results/"
save_name = save_path*content*"_L=$(L)_hmax=$(h_vals[end])_"*time_finished*".jld"

save(save_name, "result", result)
println("finished")
