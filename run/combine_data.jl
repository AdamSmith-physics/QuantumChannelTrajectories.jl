using Printf
using QuantumChannelTrajectories

#### Copy and paste from the file you ran! ###
dt = 0.25
p = 0.5
Nx = 4
Ny = 4
N = Nx*Ny
V = 1.5
b = -0.5 #2/((Nx-1)*(Ny-1))  # Magnetic field strength
num_iterations = 25
steps = 400
site_in = 1  # Site where the current is injected
drive_type = :current  # :current, :dephasing
initial_state = :random  # :checkerboard, :empty, :filled, :random, :custom
fermions = true  # Whether to use fermionic statistics
B = b*pi # Magnetic field in units of flux quantum
site_out = N  # Site where the current is extracted
###############################################

bonds = get_bonds(Nx, Ny, site_in, site_out)

K_list_accumulated = zeros(Int, steps, 9)
n_list_accumulated = zeros(Float64, steps+1, N)
currents_list_accumulated = zeros(Float64, steps+1, length(bonds))
density_correlations_accumulated = zeros(Float64, N, N)


K_avg = zeros(Int, steps, 9)
n_avg = zeros(Float64, steps+1, N)
avg_currents = zeros(Float64, steps+1, length(bonds))
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

num_processes = 40
for run_idx in 1:num_processes

    if !isfile(filename * "_run$(run_idx)" * ".h5")
        println("Skipping run $(run_idx), file does not exist with name: $(filename * "_run$(run_idx)" * ".h5")")
        continue
    end

    data = load_from_hdf5(filename * "_run$(run_idx)" * ".h5")

    global K_avg += data[:K_avg]
    global n_avg += data[:n_avg]
    global avg_currents += data[:avg_currents]
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
        :avg_currents => avg_currents ./ completed_trajectories,
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

save_to_hdf5(final_data, filename * ".h5")

# # delete files after combining
# for run_idx in 1:num_processes
#     filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))_run$(run_idx).h5"
#     rm(filename)
# end
