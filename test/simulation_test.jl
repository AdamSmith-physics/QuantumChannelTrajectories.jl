using Revise
using BenchmarkTools
using LinearAlgebra

using QuantumChannelTrajectories


Nx = 3
Ny = 4

N = Nx * Ny

num_iterations = 10

parameters::SimulationParameters = SimulationParameters(
    steps=100, 
    Nx=Nx, 
    Ny=Ny, 
    dt=0.1,
    p=0.5, 
    B=1.0,
    bonds=get_bonds(Nx, Ny, 1, 2), 
    site_in=1, 
    site_out=2, 
    drive_type="current", 
    initial_state="random")

hamiltonian = create_hamiltonian(Nx, Ny; fermions=true);
GC.gc();

ψ = random_state(Nx, Ny);

data = trajectory(hamiltonian, ψ, num_iterations, parameters)