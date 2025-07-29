
export c_op
export c_dagger_op
export n_op
export one_minus_n_op


c_local = sparse([0.0 1.0; 0.0 0.0]) # Creation operator
c_dagger_local = sparse([0.0 0.0; 1.0 0.0]) # Annihilation operator
n_local = sparse([0.0 0.0; 0.0 1.0]) # Number operator
one_minus_n_local = sparse([1.0 0.0; 0.0 0.0]) # Operator for 1 - n

function _pad_operator(op::SparseMatrixCSC{Complex{Float64}, Int}, n::Int, N::Int)
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


function pick_kraus()

end
