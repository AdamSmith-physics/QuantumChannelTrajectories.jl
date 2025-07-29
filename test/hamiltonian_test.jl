using Revise
using BenchmarkTools
using SparseArrays
using LinearAlgebra

using QuantumChannelTrajectories



@time hamiltonian = create_hamiltonian(5, 5; fermions=true)
GC.gc()
# hamiltonian
# hamiltonian = nothing

print("Hamiltonian created!")

# Wait for 20 seconds
sleep(20)