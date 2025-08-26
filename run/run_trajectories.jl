using Pkg
Pkg.activate("~/julia_projects/QuantumChannelTrajectories.jl/")
Pkg.instantiate()

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
D_list = Any[0.1, 0.25, 0.4]
# P_list = Any[0.2, 0.5]
L_list = Any[4]
V_list = Any[0.0, 2.0]
B_list = Any[0.0]
N_list = Any[500]
T_list = Any[100]
G_list = Any[false, true]


# All_input_combinations = [(d,l,v,b,n,t,g) for d in D_list, for l in L_list, for v in V_list, for b in B_list, for n in N_list, for t in T_list, for g in G_list]
All_input_combinations = [(d, l, v, b, n, t, g) for d in D_list, l in L_list, v in V_list, b in B_list, n in N_list, t in T_list, g in G_list]
run_index = (run_id % length(All_input_combinations)) +1
input_D, input_L, input_V, input_B, input_N, input_T, input_G = All_input_combinations[run_index]


dt = input_D  # Time step
p = 2*dt # Probability of hopping
Nx = input_L  # Number of sites in x-direction
Ny = input_L  # Number of sites in y-direction
V = input_V  # Interaction strength
b = input_B  # Magnetic field strength
num_iterations = input_N  # Number of iterations
steps = input_T  # Number of steps in each iteration
fermions = input_G  # Whether to use fermionic statistics


### ### ### Change these parameters as needed ### ### ###
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
trotter_evolution = true


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
    single_shot=single_shot,
    trotter_evolution = trotter_evolution
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
if trotter_evolution
    filename *= "_trotter"
end
if run_id !== nothing
    filename *= "_run$(run_id)"
end
filename *= ".h5"


if trotter_evolution
    hamiltonian = create_circuit(Nx, Ny; B=B, V=V, fermions=fermions);
else
    hamiltonian = create_hamiltonian(Nx, Ny; B=B, V=V, fermions=fermions);
end

GC.gc();

ψ = generate_initial_state(Nx, Ny; initial_state=initial_state);

run_trajectories(hamiltonian, ψ, num_iterations, fermions, parameters; eager_saving=true, filename=filename)
