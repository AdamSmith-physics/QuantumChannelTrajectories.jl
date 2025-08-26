
export create_hamiltonian
export create_circuit
export get_bonds

# BLAS.set_num_threads(1) 

"""
    create_hamiltonian(Nx::Int, Ny::Int)

    Creates a Hamiltonian for a 2D system with dimensions Nx and Ny.
"""
function create_hamiltonian(Nx::Int, Ny::Int; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)

    N::Int = Nx * Ny

    hamiltonian = spzeros(Complex{Float64}, 2^N, 2^N)

    
    row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
    for ny in 1:Ny
        local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)

        # Construct horizontal operators first
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)

        for nx in 1:Nx-1
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1))
            )
        end

        hamiltonian += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1)))
        )

    end

    

    # This does not contain magnetic field yet!
    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    local_operator = -kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
    local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
    local_operator += V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)

    for ny in 1:Ny-1
        # Construct vertical operators
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))

        for nx in 1:Nx
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx), 2^(Nx-nx))
            )
        end

        hamiltonian += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1)))
        )

    end

    row_operator = nothing
    local_operator = nothing

    return hamiltonian

end


function create_circuit(Nx::Int, Ny::Int; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)

    N::Int = Nx * Ny

    layer_1 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_2 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_3 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_4 = spzeros(Complex{Float64}, 2^N, 2^N)

    row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
    for ny in 1:Ny
        local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)

        # Constructing first layer of horizontal operators
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
        for nx in 1:2:Nx-1
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1))
            )
        end
        layer_1 += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1)))
            )        
                
        # Constructing second layer of horizontal operators
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)    
        for nx in 2:2:Nx-1
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1))
            )
        end

        layer_2 += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1)))
        )

    end    

    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)

    local_operator = V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
    local_operator -= kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
    local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
    
    # creating first layer of vertical operators    
    for ny in 1:2:Ny-1   
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        for nx in 1:Nx
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx), 2^(Nx-nx))
            )
        end

        layer_3 += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1)))
        )
    end
    
    # creating second layer of vertical operators
    for ny in 2:2:Ny-1 
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        
        for nx in 1:Nx
            row_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx), 2^(Nx-nx))
                )
        end
            
        layer_4 += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            row_operator,
            sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1)))
            )            
    end
            
            row_operator = nothing
            local_operator = nothing
            
    return layer_1, layer_2, layer_3, layer_4

end


function get_bonds(Nx::Int, Ny::Int, site_in::Int, site_out::Int)::Vector{Tuple{Int,Int}}
    bonds = []

    # Horizontal bonds
    for ny in 1:Ny
        for nx in 1:Nx-1
            n1 = nx + (ny-1)*Nx
            n2 = (nx+1) + (ny-1)*Nx

            if !( n1 ∈ [site_in, site_out] ) && !( n2 ∈ [site_in, site_out] )
                push!(bonds, (n1, n2))
            end
        end
    end

    # Vertical bonds
    for ny in 1:Ny-1
        for nx in 1:Nx
            n1 = nx + (ny-1)*Nx
            n2 = nx + ny*Nx

            if !( n1 ∈ [site_in, site_out] ) && !( n2 ∈ [site_in, site_out] )
                push!(bonds, (n1, n2))
            end
        end
    end

    return bonds
end
