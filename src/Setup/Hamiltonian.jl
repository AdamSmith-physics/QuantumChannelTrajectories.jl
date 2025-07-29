
export create_hamiltonian

# BLAS.set_num_threads(1) 

PauliX = sparse([0.0 1.0; 1.0 0.0])
PauliY = sparse([0.0 -im; im 0.0])
PauliZ = sparse([1.0 0.0; 0.0 -1.0])
density_operator = sparse([0.0 0.0; 0.0 1.0])
PauliOperators = [PauliX, PauliY, PauliZ]


"""
    create_hamiltonian(Nx::Int, Ny::Int)

    Creates a Hamiltonian for a 2D system with dimensions Nx and Ny.
"""
function create_hamiltonian(Nx::Int, Ny::Int; V::Float64 = 0.0, fermions::Bool = false)::SparseMatrixCSC{Complex{Float64}, Int}

    N::Int = Nx * Ny

    hamiltonian = spzeros(Complex{Float64}, 2^N, 2^N)

    local_operator = -kron(PauliX,PauliX) - kron(PauliY,PauliY) + V*kron(PauliZ,PauliZ)
    for ny in 1:Ny
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

    row_operator = spzeros(Complex{Float64}, 1, 1)
    # GC.gc()

    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    local_operator = -kron(PauliX, fill(fill_operator,Nx-1)... ,PauliX) / 2
    local_operator -= kron(PauliY, fill(fill_operator,Nx-1)..., PauliY) / 2
    local_operator += V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
    for ny in 1:Ny-1
        # Construct vertical operators
        col_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))

        for nx in 1:Nx
            col_operator += kron(
                sparse(I, 2^(nx-1), 2^(nx-1)),
                local_operator,
                sparse(I, 2^(Nx-nx), 2^(Nx-nx))
            )
        end

        hamiltonian += kron(
            sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),
            col_operator,
            sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1)))
        )

    end

    col_operator = spzeros(Complex{Float64}, 1, 1)
    # GC.gc()

    local_operator = spzeros(Complex{Float64}, 1, 1)
    # GC.gc()

    return hamiltonian

end


function get_bonds(Nx::Int, Ny::Int, site_in::Int, site_out::Int)::Vector{Tuple{Int, Int}}
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
