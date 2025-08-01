
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
    (x1, y1) = ((n1-1) % Nx + 1, (n1-1) ÷ Nx + 1)
    (x2, y2) = ((n2-1) % Nx + 1, (n2-1) ÷ Nx + 1)

    res = 0.0

    if x1 == x2  # Vertical bond

        res += -im * ψ' * (_c_dag_c(n1,n2,N,fermions) * ψ)
        res += im * ψ' * (_c_dag_c(n2,n1,N,fermions) * ψ)
    
    elseif y1 == y2  # Horizontal bond

        res += -im * exp(im*B*y1) * ψ' * (_c_dag_c(n1,n2,N,fermions)  * ψ)
        res += im * exp(-im*B*y1) * ψ' * (_c_dag_c(n2,n1,N,fermions) * ψ)

    end

    return real(res)

end


function _c_dag_c(n1::Int, n2::Int, N::Int, fermions::Bool)
    # More efficient calculation of <c_dagger(n1) c(n2)>

    if n1 == n2
        error("Shouldn't be called with n1 == n2. This is for density!")
    end

    op1 = Sigma_plus
    op2 = Sigma_minus
    if n2 < n1
        op1 = Sigma_minus
        op2 = Sigma_plus
    end

    n_min = min(n1, n2)
    n_max = max(n1, n2)

    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    return kron(
        [[sparse(I, 2^(n_min-1), 2^(n_min-1)),
        op1];
        fill(fill_operator, n_max-n_min-1);
        [op2,
        sparse(I, 2^(N - n_max), 2^(N - n_max))]]...
    )
end
