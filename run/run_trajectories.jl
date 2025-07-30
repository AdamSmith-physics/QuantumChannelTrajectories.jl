using Revise
using LinearAlgebra

using QuantumChannelTrajectories

BLAS.set_num_threads(1) 

run_idx = nothing
if length(ARGS) > 0
    run_idx = parse(Int, ARGS[1])
end

print(run_idx)
typeof(run_idx)

dt = 0.2
p = 0.2
Nx = 4
Ny = 4
N = Nx*Ny
V = 0.0
b = 0.0 #2/((Nx-1)*(Ny-1))  # Magnetic field strength
B = b*pi # Magnetic field in units of flux quantum
num_iterations = 100
steps = 25
site_in = 1  # Site where the current is injected
site_out = N  # Site where the current is extracted
drive_type = :current  # :current, :dephasing
initial_state = :random  # :checkerboard, :empty, :filled, :random, :custom
fermions = true  # Whether to use fermionic statistics

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

hamiltonian = create_hamiltonian(Nx, Ny; fermions=true);
GC.gc();

ψ = generate_initial_state(Nx, Ny; initial_state=initial_state);

data = run_trajectories(hamiltonian, ψ, num_iterations, fermions, parameters)

using Printf
filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))"
if run_idx !== nothing
    filename *= "_run$(run_idx)"
end
filename *= ".h5"
save_to_hdf5(data, filename)