using Revise
using BenchmarkTools
using LinearAlgebra

using QuantumChannelTrajectories

BLAS.set_num_threads(1) 

Nx = 3
Ny = 3

N = Nx * Ny

num_iterations = 100

parameters = SimulationParameters(
    steps=200, 
    Nx=Nx, 
    Ny=Ny, 
    dt=0.1,
    p=0.5, 
    B=1.0,
    bonds=get_bonds(Nx, Ny, 1, 2), 
    site_in=1, 
    site_out=2, 
    drive_type=:current, 
    initial_state=:random)


hamiltonian = create_hamiltonian(Nx, Ny; fermions=true);
GC.gc();

ψ = random_state(Nx, Ny);

data = run_trajectories(hamiltonian, ψ, num_iterations, parameters)

save_to_hdf5(data, "simulation_results.h5")

# Load the data back
loaded_data = load_from_hdf5("simulation_results.h5")

test = load_key_from_hdf5("simulation_results.h5", :parameters)

load_parameters = loaded_data[:parameters]

data[:parameters][:bonds]
loaded_data[:parameters][:bonds]

data == loaded_data

data
loaded_data

data[:K_list] == loaded_data[:K_list]
data == loaded_data