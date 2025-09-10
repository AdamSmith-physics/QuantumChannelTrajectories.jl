using Printf
using QuantumChannelTrajectories


dt = 0.25
p = 0.5
Nx = 4
Ny = 4
N = Nx*Ny
V = 4.0
b = 0.0 #2/((Nx-1)*(Ny-1))  # Magnetic field strength
num_iterations = 100
steps = 100
site_in = 1  # Site where the current is injected
drive_type = :current  # :current, :dephasing
initial_state = :custom  # :checkerboard, :empty, :filled, :random, :custom
fermions = false  # Whether to use fermionic statistics
B = b*pi # Magnetic field in units of flux quantum
site_out = N  # Site where the current is extracted
# Optional parameters
even_parity = false  # Whether to enforce even parity
pinned_corners = true  # Whether to pin the corners
single_shot = false  # Whether to perform single shot measurements
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

num_processes = 15
for run_idx in 1:num_processes

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

    if run_idx == 1
        global t_list = data[:t_list]
        global parameters = data[:params]
    end

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

save_to_hdf5(final_data, filename * ".h5")

# # delete files after combining
# for run_idx in 1:num_processes
#     filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))_run$(run_idx).h5"
#     rm(filename)
# end
