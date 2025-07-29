using Revise
using BenchmarkTools
using SparseArrays
using LinearAlgebra
using KrylovKit

using QuantumChannelTrajectories

BLAS.set_num_threads(1) 

Nx = 5
Ny = 5


hamiltonian = create_hamiltonian(Nx, Ny; fermions=true);
GC.gc();
# hamiltonian
# hamiltonian = nothing

# create a random product state
ψ = product_state(rand(0:1, Nx * Ny));

@benchmark hamiltonian * ψ

@benchmark ψ = exponentiate($hamiltonian, im*1.01, $ψ; tol=1e-14)

@time ψ, info = exponentiate(hamiltonian, -im*0.1, ψ; tol=1e-15, ishermitian=true, eager=true)   
println(ψ)
print("Hamiltonian created!")

# Wait for 20 seconds
sleep(20)