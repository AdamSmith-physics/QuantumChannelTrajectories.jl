using Pkg
Pkg.activate("~/julia_projects/QuantumChannelTrajectories.jl/")
Pkg.instantiate()

using Printf
include("../src/QuantumChannelTrajectories.jl")
using .QuantumChannelTrajectories


run_id = 1  
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


B = b*pi # Magnetic field in units of flux quantum
N = Nx*Ny
site_in = 1  # Site where the current is injected
site_out = N  # Site where the current is extracted
drive_type = :current  # :current, :dephasing
initial_state = :random  # :checkerboard, :empty, :filled, :random, :custom
# Optional parameters
even_parity = false  # Whether to enforce even parity
pinned_corners = true  # Whether to pin the corners
single_shot = false  # Whether to perform single shot measurements
trotter_evolution = true
###############################################

bonds = get_bonds(Nx, Ny, site_in, site_out)


K_avg = zeros(Int, steps, 9)
n_avg = zeros(Float64, steps+1, N)
n_sq_avg = zeros(Float64, steps+1, N)
avg_currents = zeros(Float64, steps+1, length(bonds))
currents_sq_avg = zeros(Float64, steps+1, length(bonds))
avg_dd_correlations = zeros(Float64, N, N)
completed_trajectories = 0
t_list = nothing
parameters = nothing

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

num_processes = 50
off_set = length(All_input_combinations)
for run_idx in run_id-1:off_set:off_set*num_processes


    if !isfile(filename * "_run$(run_idx)" * ".h5")
        println("Skipping run $(run_idx), file does not exist with name: $(filename * "_run$(run_idx)" * ".h5")")
        continue
    end

    data = load_from_hdf5(filename * "_run$(run_idx)" * ".h5")

    global K_avg += data[:K_avg]
    global n_avg += data[:n_avg]
    global n_sq_avg += data[:n_sq_avg]
    global avg_currents += data[:avg_currents]
    global currents_sq_avg += data[:currents_sq_avg]
    global avg_dd_correlations += data[:avg_dd_correlations]
    global completed_trajectories += data[:completed_trajectories]

    # if run_idx == 1
    global t_list = data[:t_list]
    global parameters = data[:params]
    # end

end

final_data = Dict(
        :K_avg => K_avg ./ completed_trajectories,
        :n_avg => n_avg ./ completed_trajectories,
        :n_sq_avg => n_sq_avg ./ completed_trajectories,
        :avg_currents => avg_currents ./ completed_trajectories,
        :currents_sq_avg => currents_sq_avg ./ completed_trajectories,
        :avg_dd_correlations => avg_dd_correlations ./ completed_trajectories,
        :t_list => t_list,
        :params => parameters
    )


filename = ""
if fermions
    filename = "data/fermions_"
else
    filename = "data/bosons_"
end
filename *= "$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_V$(V)_steps$(steps)_trajectories$(completed_trajectories)_$(string(drive_type))_$(string(initial_state))"
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
save_to_hdf5(final_data, filename * ".h5")

# # delete files after combining
# for run_idx in 1:num_processes
#     filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))_run$(run_idx).h5"
#     rm(filename)
# end
