
export product_state
export checkerboard_state
export empty_state
export filled_state
export random_state
export generate_initial_state


# include("/Users/arash/Julia Projects/Quantum_Trajectory_Adam/QuantumChannelTrajectories.jl/src/Other/Measurement.jl")
# using .Measurement.jl
# using Measurement.jl
# using .QuantumChannelTrajectories.Measurement


"""
    product_state(occupation_list::Vector{Int})

Create a product state from a list of occupations.

# Arguments
- `occupation_list`: A vector of integers representing the occupation numbers for each mode.

# Returns
A vector representing the product state.
"""

function single_shot(value::Float64, min::Float64, max::Float64)::Float64
    sample = (max-min)*rand() + min
    if sample < value
        return max
    else
        return min
    end
end

function product_state(occupation_list::Vector{Int}, Nx::Int, Ny::Int)
    
    # Check if the occupation list is the right length
    if length(occupation_list) != Nx * Ny
        throw(ArgumentError("Occupation list must have length Nx * Ny"))
    end

    ψ::Array{Complex{Float64},1} = [1.0]

    # loop over enumeration of occupation_list
    for occupation in occupation_list
        # create a local state for each occupation
        new_state = (occupation == 0) ? [1.0, 0.0] : [0.0, 1.0]

        ψ = kron(ψ, new_state)
    end

    return ψ
end


function checkerboard_state(Nx::Int, Ny::Int)
    occupation_list = zeros(Int, Nx * Ny)

    # Fill the occupation list with a checkerboard pattern
    for ny in 1:Ny
        for nx in 1:Nx
            if (nx + ny) % 2 == 1
                occupation_list[nx + (ny-1)*Nx] = 1
            end
        end
    end

    return product_state(occupation_list, Nx, Ny)
end


function empty_state(Nx::Int, Ny::Int)
    occupation_list = zeros(Int, Nx * Ny)
    return product_state(occupation_list, Nx, Ny)
end


function filled_state(Nx::Int, Ny::Int)
    occupation_list = ones(Int, Nx * Ny)
    return product_state(occupation_list, Nx, Ny)
end


function random_state(Nx::Int, Ny::Int; even_parity::Bool = false, pinned_corners::Bool = false, site_in::Int = 1, site_out::Int = Nx * Ny, initialization::Any=nothing)
    
    occupation_list = rand(0:1, Nx * Ny)
    
    if initialization !== nothing
        occupation_list = [ Int.(single_shot(n_ini, 0.0, 1.0)) for n_ini in initialization]
    end


    if pinned_corners
        occupation_list[site_in] = 1
        occupation_list[site_out] = 0
    end

    # fix if even parity is required
    if even_parity
        # create list of sites not including site_in and site_out
        sites = collect(1:Nx*Ny)
        sites = deleteat!(sites, [site_in, site_out])
        # select a site from the list at random
        site = rand(sites)
        occupation_list[site] = (sum(occupation_list) - occupation_list[site]) % 2
    end

    return product_state(occupation_list, Nx, Ny)
end


function generate_initial_state(Nx::Int, Ny::Int; initial_state::Symbol = :random, initialization::Any=nothing)
    if initial_state == :checkerboard
        return checkerboard_state(Nx, Ny)
    elseif initial_state == :empty
        return empty_state(Nx, Ny)
    elseif initial_state == :filled
        return filled_state(Nx, Ny)
    elseif initial_state == :random
        return random_state(Nx, Ny; initialization=initialization)  # optional parameters not needed since they will be overwritten each trajectory
    else
        throw(ArgumentError("Unknown initial state: $initial_state"))
    end
end