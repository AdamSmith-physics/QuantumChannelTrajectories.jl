module QuantumChannelTrajectories

using LinearAlgebra
using SparseArrays
using KrylovKit

include("Setup/Hamiltonian.jl")
include("Setup/InitialState.jl")

include("Simulation/KrausOperators.jl")
include("Simulation/Trajectories.jl")

end # module QuantumChannelTrajectories
