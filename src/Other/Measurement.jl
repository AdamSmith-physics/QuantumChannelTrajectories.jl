
export n_expectation
export current_expectation
export density_correlations


function n_expectation(ψ::Vector{Complex{Float64}}, site::Int, N::Int; single_shot::Bool = false)
    # Calculate the expectation value of the number operator at a given site
    n_op_matrix = n_op(site, N)
    n_exp = real(ψ' * n_op_matrix * ψ)

    if !single_shot
        return n_exp
    else
        return _single_shot(n_exp, 0.0, 1.0)
    end
end


function current_expectation(ψ::Vector{Complex{Float64}}, B::Float64, bond::Tuple{Int,Int}, Nx::Int, Ny::Int, fermions::Bool; single_shot::Bool = false)
    # Calculate the expectation value of the current operator for a given bond
    
    N = Nx * Ny
    
    n1, n2 = bond
    (x1, y1) = ((n1-1) % Nx + 1, (n1-1) ÷ Nx + 1)
    (x2, y2) = ((n2-1) % Nx + 1, (n2-1) ÷ Nx + 1)

    res = 0.0
    term1 = 0.0
    term2 = 0.0

    if x1 == x2  # Vertical bond

        term1 = real(-im * ψ' * (_c_dag_c(n1,n2,N,fermions) * ψ))
        term2 = real(im * ψ' * (_c_dag_c(n2,n1,N,fermions) * ψ))
    
    elseif y1 == y2  # Horizontal bond

        term1 = real(-im * exp(im*B*y1) * ψ' * (_c_dag_c(n1,n2,N,fermions)  * ψ))
        term2 = real(im * exp(-im*B*y1) * ψ' * (_c_dag_c(n2,n1,N,fermions) * ψ))

    end

    if single_shot
        # Need to sample each term independently
        term1 = _single_shot(term1, -1.0, 1.0)
        term2 = _single_shot(term2, -1.0, 1.0)
    end

    res = term1 + term2

    return res

end


function density_correlations(ψ::Vector{Complex{Float64}}, Nx::Int, Ny::Int)
    # Doesn't currently support single shot!

    N = Nx * Ny

    density_correlations = zeros(Float64, N, N)

    for n1 in 1:N
        for n2 in n1+1:N
            # Calculate <n(n1) n(n2)> - <n(n1)> <n(n2)>
            # This is the connected density-density correlation function

            density_correlations[n1, n2] = real(ψ' * _n_n(n1, n2, N) * ψ) - n_expectation(ψ, n1, N) * n_expectation(ψ, n2, N)
        end
    end

    density_correlations += density_correlations'  # Make it symmetric

    for n in 1:N
        density_correlations[n, n] = n_expectation(ψ, n, N) - n_expectation(ψ, n, N)^2
    end

    return density_correlations
end


function _single_shot(value::Float64, min::Float64, max::Float64)::Float64

    sample = (max-min)*rand() + min

    if sample < value
        return max
    else
        return min
    end
end


function _n_n(n1::Int, n2::Int, N::Int)
    # More efficient calculation of <n(n1) n(n2)>
    
    n_min = min(n1, n2)
    n_max = max(n1, n2)

    return kron(
        sparse(I, 2^(n_min-1), 2^(n_min-1)),
        density_operator,
        sparse(I, 2^(n_max-n_min-1), 2^(n_max-n_min-1)),
        density_operator,
        sparse(I, 2^(N - n_max), 2^(N - n_max))
    )

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
