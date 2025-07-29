
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

    if K == 0
        # No Kraus operator applied
        return ψ
    end


    This needs completing! Getting a bit of a mess!
    if drive_type == "current"
        # Apply current Kraus operator
        if K == 1
            # Apply inflow Kraus operator
            ψ = c_dagger_op(site, N) * ψ
        elseif K == 2
            # Apply outflow Kraus operator
            ψ = n_op(site, N) * ψ
        end
   
    elseif drive_type == "dephasing"
        # Apply dephasing Kraus operator
        if K == 1
            # Apply 1-n
            ψ = one_minus_n_op(site, N) * ψ
        elseif K == 2
            # Apply n
            ψ = n_op(site, N) * ψ
        end
    end








    return normalize!(ψ)
end
