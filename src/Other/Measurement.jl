
export n_expectation

function n_expectation(ψ::Vector{Complex{Float64}}, site::Int, N::Int)
    # Calculate the expectation value of the number operator at a given site
    n_op_matrix = n_op(site, N)
    return real(ψ' * n_op_matrix * ψ)
end


function current_expectation(ψ::Vector{Complex{Float64}}, B::Float64, bond::Tuple{Int,Int}, Nx::Int, Ny::Int, fermions::Bool)
    # Calculate the expectation value of the current operator for a given bond
    
    N = Nx * Ny
    
    n1, n2 = bond
    (x1, y1) = ((n1-1) % Nx + 1, n1 ÷ Nx + 1)
    (x2, y2) = ((n2-1) % Nx + 1, n2 ÷ Nx + 1)
    
    res = 0.0

    if x1 == x2  # Vertical bond

        res += -im * ψ' * (c_dagger_op(n1, N, fermions) * (c_op(n2, N, fermions) * ψ))
        res += im * ψ' * (c_dagger_op(n2, N, fermions) * (c_op(n1, N, fermions) * ψ))
    
    elseif y1 == y2  # Horizontal bond

        res += -im * exp(im*B*y1) * ψ' * (c_dagger_op(n1, N, fermions) * (c_op(n2, N, fermions) * ψ))
        res += im * exp(-im*B*y1) * ψ' * (c_dagger_op(n2, N, fermions) * (c_op(n1, N, fermions) * ψ))

    end

    return real(res)

end