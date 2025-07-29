using Revise
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using ExponentialAction
using ExponentialUtilities
using KrylovKit

using QuantumChannelTrajectories

# Create a 1D array of 1000 complex elements
arr = rand(Complex{Float64}, 2^25);

# Compute norm of the array
@benchmark norm($arr)


M = rand(Complex{Float64}, 1000, 1000);

M * M
M

# Create empty sparse matrix
S = spzeros(Complex{Float64}, 1000, 1000)

X = sparse([0.0 1.0; 1.0 0.0])
Y = sparse([0.0 -im; im 0.0])
Z = sparse([1.0 0.0; 0.0 -1.0])

operator_list = [X, Y, Z]
# Use kron to create a tensor product of the operators
M = kron(operator_list...)


M = kron(X, Z, Z, X)

test = create_hamiltonian(2, 2)

# List containing Z, 10 times
Z_list = fill(Z, 25)
M = kron(Z_list...)
M += kron(fill(X, 25)...)
M += kron(fill(Y, 25)...)

BLAS.set_num_threads(1)


@benchmark ExponentialAction.expv(im*1.01, $M, $arr, tol=1e-15)
@benchmark ExponentialUtilities.expv(im*1.01, $M, $arr, tol=1e-15)
@benchmark KrylovKit.exponentiate($M, im*1.01, $arr; tol=1e-15)
