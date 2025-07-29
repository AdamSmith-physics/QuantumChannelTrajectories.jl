
export product_state
export checkerboard_state
export empty_state
export filled_state
export random_state

"""
    product_state(occupation_list::Vector{Int})

Create a product state from a list of occupations.

# Arguments
- `occupation_list`: A vector of integers representing the occupation numbers for each mode.

# Returns
A vector representing the product state.
"""
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
    for nx in 1:Nx
        for ny in 1:Ny
            if (i + j) % 2 == 0
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


function random_state(Nx::Int, Ny::Int; even_parity::Bool = false)
    occupation_list = rand(0:1, Nx * Ny)

    # fix if even parity is required
    if even_parity
        occupation_list[1] = sum(occupation_list[2:end]) % 2
    end

    return product_state(occupation_list, Nx, Ny)
end