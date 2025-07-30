
export c_op
export c_dagger_op
export n_op
export one_minus_n_op

export pick_kraus


c_local = sparse([0.0 1.0; 0.0 0.0]) # Creation operator
c_dagger_local = sparse([0.0 0.0; 1.0 0.0]) # Annihilation operator
n_local = sparse([0.0 0.0; 0.0 1.0]) # Number operator
one_minus_n_local = sparse([1.0 0.0; 0.0 0.0]) # Operator for 1 - n

function _pad_operator(op::SparseMatrixCSC, n::Int, N::Int)
    # Pad the operator to the full size of the system
    return kron(sparse(I, 2^(n-1), 2^(n-1)), op, sparse(I, 2^(N - n), 2^(N - n)))
end


function c_op(n::Int, N::Int)
    # Create a local creation operator for site n
    return _pad_operator(c_local, n, N)
end


function c_dagger_op(n::Int, N::Int)
    # Create a local annihilation operator for site n
    return _pad_operator(c_dagger_local, n, N)
end


function n_op(n::Int, N::Int)
    # Create a local number operator for site n
    return _pad_operator(n_local, n, N)
end


function one_minus_n_op(n::Int, N::Int)
    # Create a local operator for 1 - n for site n
    return _pad_operator(one_minus_n_local, n, N)
end


function pick_kraus(p::Float64, ψ::Vector{Complex{Float64}}, site::Int, N::Int)
    # Pick index of the Kraus operator from 0 to 2
    
    n_exp = n_expectation(ψ, site, N)

    probabilities = cumsum([1-p, p*(1-n_exp), p*n_exp])
 
    return sum(rand() .> probabilities)

end


function apply_kraus!(ψ::Vector{Complex{Float64}}, K::Int, site::Int, N::Int; type=:inflow, drive_type="current")

    # This top part is same for current/dephasing and inflow/outflow
    if K == 0
        # No Kraus operator applied
        return nothing
    elseif K == 1
        ψ = one_minus_n_op(site, N) * ψ
    elseif K == 2
        ψ = n_op(site, N) * ψ
    else
        error("Invalid Kraus operator index: $K")
    end

    if drive_type == :current && type == :inflow && K == 1
        # Apply inflow Kraus operator
        ψ = c_dagger_op(site, N) * ψ
    elseif drive_type == :current && type == :outflow && K == 2
        # Apply outflow Kraus operator
        ψ = c_op(site, N) * ψ
    end

    normalize!(ψ)

    return nothing
end
