using Revise
using LinearAlgebra
using Printf

using QuantumChannelTrajectories

BLAS.set_num_threads(1) 

run_id = 1  # Default run_id = 1. The results will be summed, not averaged. Use combine_data.jl to average results across multiple runs.
if length(ARGS) > 0
    run_id = parse(Int, ARGS[1])
end

dt = 0.25
p = 0.25
Nx = 4
Ny = 4
N = Nx*Ny
V = 0.0
b = 0.0 #2/((Nx-1)*(Ny-1))  # Magnetic field strength
num_iterations = 30
steps = 50
site_in = 1  # Site where the current is injected
drive_type = :current  # :current, :dephasing
initial_state = :random  # :checkerboard, :empty, :filled, :random, :custom
fermions = true  # Whether to use fermionic statistics
B = b*pi # Magnetic field in units of flux quantum
site_out = N  # Site where the current is extracted


parameters = SimulationParameters(
    steps=steps, 
    Nx=Nx, 
    Ny=Ny, 
    dt=dt,
    p=p, 
    B=B,
    bonds=get_bonds(Nx, Ny, site_in, site_out), 
    site_in=site_in, 
    site_out=site_out, 
    drive_type=drive_type, 
    initial_state=initial_state
    )


filename = ""
if fermions
    filename = "data/fermions_"
else
    filename = "data/bosons_"
end

filename *= "$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_V$(V)_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))"
if run_id !== nothing
    filename *= "_run$(run_id)"
end
filename *= ".h5"


hamiltonian = create_hamiltonian(Nx, Ny; B=B, V=V, fermions=fermions);
GC.gc();

ψ = generate_initial_state(Nx, Ny; initial_state=initial_state);

run_trajectories(hamiltonian, ψ, num_iterations, fermions, parameters; eager_saving=true, filename=filename)
