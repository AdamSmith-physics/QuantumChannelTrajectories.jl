using Revise
using LinearAlgebra
using Printf

using QuantumChannelTrajectories

BLAS.set_num_threads(1) 

run_id = 1  # Default run_id = 1. The results will be summed, not averaged. Use combine_data.jl to average results across multiple runs.
if length(ARGS) > 0
    run_id = parse(Int, ARGS[1])
end

# Parameters set at runtime
dt = parse(Float64, ARGS[2])  # Time step
p = parse(Float64, ARGS[3])  # Probability of hopping
Nx = parse(Int, ARGS[4])  # Number of sites in x-direction
Ny = parse(Int, ARGS[5])  # Number of sites in y-direction
V = parse(Float64, ARGS[6])  # Interaction strength
b = parse(Float64, ARGS[7])  # Magnetic field strength
num_iterations = parse(Int, ARGS[8])  # Number of iterations
steps = parse(Int, ARGS[9])  # Number of steps in each iteration
fermions = parse(Bool, ARGS[10])  # Whether to use fermionic statistics

# Change these parameters as needed
N = Nx*Ny
site_in = 1  # Site where the current is injected
drive_type = :current  # :current, :dephasing
initial_state = :random  # :checkerboard, :empty, :filled, :random, :custom
B = b*pi # Magnetic field in units of flux quantum
site_out = N  # Site where the current is extracted

# Optional parameters
even_parity = false  # Whether to enforce even parity
pinned_corners = true  # Whether to pin the corners
single_shot = true


println("\nRunning with parameters:")
# print all parameters separated by newline
println("dt: $dt \n",
        "p: $p \n",
        "Nx: $Nx \n",
        "Ny: $Ny \n",
        "V: $V \n",
        "b: $b \n",
        "num_iterations: $num_iterations \n",
        "steps: $steps \n",
        "fermions: $fermions \n",
        "site_in: $site_in \n",
        "drive_type: $drive_type \n",
        "initial_state: $initial_state \n",
        "B: $B \n",
        "site_out: $site_out \n")


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
    initial_state=initial_state,
    even_parity=even_parity,
    pinned_corners=pinned_corners,
    single_shot=single_shot
    )



filename = ""
if fermions
    filename = "data/fermions_"
else
    filename = "data/bosons_"
end

filename *= "$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_V$(V)_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))"
if even_parity
    filename *= "_even_parity"
end
if pinned_corners
    filename *= "_pinned_corners"
end
if single_shot
    filename *= "_single_shot"
end
if run_id !== nothing
    filename *= "_run$(run_id)"
end
filename *= ".h5"


hamiltonian = create_hamiltonian(Nx, Ny; B=B, V=V, fermions=fermions);
GC.gc();

ψ = generate_initial_state(Nx, Ny; initial_state=initial_state);

run_trajectories(hamiltonian, ψ, num_iterations, fermions, parameters; eager_saving=true, filename=filename)
