export create_circuit
export create_circuit_Adam
export default_circuit_order


""" Checks if the tuples in each layer are disjoint pairs of site indices. I.e., no site index appears more than once in a layer.
    Also checks that the site indices are within the valid range [1, Nx*Ny], and that the tuples are nearest-neighbours.
    
    Returns true if the order is valid, false otherwise.
"""
function _check_trotter_order(Nx::Int, Ny::Int, order)::Tuple{Bool,String}
    
    for layer in order
        seen_indices = Set{Int}()
        for (i, j) in layer
            if i in seen_indices || j in seen_indices
                return false, "Overlapping indices in layer detected."
            end
            push!(seen_indices, i)
            push!(seen_indices, j)
        end
    end

    N = Nx * Ny
    for layer in order
        for (i, j) in layer
            if i < 1 || i > N || j < 1 || j > N
                return false, "Indices ($i, $j) out of bounds."  # Indices out of bounds
            end
            if abs(i - j) != 1 && abs(i - j) != Nx
                return false, "Indices ($i, $j) are not nearest-neighbours."  # Not nearest neighbours
            end
        end
    end

    return true, "Trotter order is valid."

end


"""
    create_circuit(Nx::Int, Ny::Int, Order)

    Creates a Trotter circuit for a 2D system with dimensions Nx and Ny,
    following the specified order of two-site interactions.

    Parameters:
    - `Nx::Int`: Number of sites in the x-direction.
    - `Ny::Int`: Number of sites in the y-direction.
    - `Order::Vector{Vector{Tuple{Int, Int}}}`: A vector of layers, where each layer is a vector of tuples representing pairs of site indices to interact.
    - `B::Float64`: Magnetic field strength (default is 0.0).
    - `V::Float64`: Interaction strength (default is 0.0).
    - `fermions::Bool`: Whether to use fermionic statistics (default is false).
"""
function create_circuit(Nx::Int, Ny::Int, order; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false, interaction_type::Symbol = :density)
    
    # Check interaction type is valid, options are :density or :ZZ
    if !(interaction_type in [:density, :ZZ])
        error("Invalid interaction type: $interaction_type. Valid options are :density or :ZZ.")
    end

    trotter_valid, message = _check_trotter_order(Nx, Ny, order)
    if !trotter_valid
        error("Invalid Trotter order: $message")
    end

    N::Int = Nx * Ny

    num_layers = length(order)
    circuit_layers = [spzeros(Complex{Float64}, 2^N, 2^N) for _ in 1:num_layers]  # Create a Hamiltonian for each layer
    
    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)  # Operator to fill the JW string with (fermions vs bosons)
    
    for layer in 1:num_layers
        num_indices = length(order[layer])
        for indices in 1:num_indices
            x1, x2 = order[layer][indices]
            nx::Int = mod(x1-1,Nx)+1
            ny::Int = div(x1-1,Nx)+1
            
            if abs(x1-x2)==1

                local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus)
                if interaction_type == :density
                    local_operator += V*kron(density_operator,density_operator)
                elseif interaction_type == :ZZ
                    local_operator += V/4*kron(PauliZ,PauliZ)  # V/4 because Sigma_z = 2n - I
                end

                row_operator = kron( 
                    sparse(I, 2^(nx-1), 2^(nx-1)), 
                    local_operator, 
                    sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)) 
                )

                circuit_layers[layer] += kron( 
                    sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), 
                    row_operator, 
                    sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1)))
                )        
            
            elseif  abs(x1-x2)==Nx

                local_operator = -kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
                local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
                if interaction_type == :density
                    local_operator += V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
                elseif interaction_type == :ZZ
                    local_operator += V*kron(PauliZ, sparse(I,2^(Nx-1),2^(Nx-1)), PauliZ)  # V/4 because Sigma_z = 2n - I
                end

                row_operator = kron( 
                    sparse(I, 2^(nx-1), 2^(nx-1)), 
                    local_operator, 
                    sparse(I, 2^(Nx-nx), 2^(Nx-nx)) 
                )

                circuit_layers[layer] += kron(
                    sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), 
                    row_operator, 
                    sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))) 
                ) 
                          
            else
                error("Indices ($x1, $x2) are not nearest-neighbours. This wasn't caught in the initial check.")
            end
        end
    end

    # manual cleanup
    row_operator = nothing
    local_operator = nothing

    return circuit_layers
end


function default_circuit_order(Nx::Int, Ny::Int; staggered::Bool, alternating::Bool)

    layers=[[] for _ in 1:4]

    index_order = alternating ? [1,2,3,4] : [1,3,2,4]

    # Horizontal bonds
    for ny in 1:Ny
        # shift 0 starts on site 1, shift 1 starts on site 2
        for shift in [0,1]
            start_col = (staggered*(ny-1)+shift)%2 + 1  # Determine starting column based on staggered and shift

            for nx in start_col:2:Nx-1
                index = nx + (ny-1)*Nx
                
                push!(layers[index_order[1+2*shift]], (index, index+1)) # Horizontal bonds
            end
        end
    end

    # Vertical bonds
    for nx in 1:Ny
        # shift 0 starts on site 1, shift 1 starts on site 2
        for shift in [0,1]
            start_col = (staggered*(nx-1)+shift)%2 + 1  # Determine starting column based on staggered and shift

            for ny in start_col:2:Ny-1
                index = nx + (ny-1)*Nx
                
                push!(layers[index_order[2+2*shift]], (index, index+Nx)) # Vertical bonds
            end
        end
    end


    return layers

end



