using Revise

using QuantumChannelTrajectories

if length(ARGS) > 0
    run = parse(Int, ARGS[1])
end

print(run)

dt = 0.2
p = 0.3
Nx = 3
Ny = 3
N = Nx*Ny
V = 0.0
b = 0.0 #2/((Nx-1)*(Ny-1))  # Magnetic field strength
B = b*pi # Magnetic field in units of flux quantum
num_iterations = 10
steps = 50
site_in = 1  # Site where the current is injected
site_out = N  # Site where the current is extracted
drive_type = :current  # "current", "dephasing"
initial_state = :random  # "checkerboard", "empty", "random", "custom"
fermions = true  # Whether to use fermionic statistics

parameters = SimulationParameters(
    steps=steps, 
    Nx=Nx, 
    Ny=Ny, 
    dt=dt,
    p=p, 
    B=B,
    bonds=get_bonds(Nx, Ny, 1, 2), 
    site_in=site_in, 
    site_out=site_out, 
    drive_type=drive_type, 
    initial_state=initial_state
    )

hamiltonian = create_hamiltonian(Nx, Ny; fermions=true);
GC.gc();

ψ = random_state(Nx, Ny);

data = run_trajectories(hamiltonian, ψ, num_iterations, fermions, parameters)

using Printf
filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_B$(@sprintf("%.2f", B))_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))_$(run).h5"
save_to_hdf5(data, filename)