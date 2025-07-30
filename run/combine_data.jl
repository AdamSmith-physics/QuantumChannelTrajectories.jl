using Printf
using QuantumChannelTrajectories

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

bonds = get_bonds(Nx, Ny, site_in, site_out)

K_list_accumulated = zeros(Int, steps, 9)
n_list_accumulated = zeros(Float64, steps+1, N)
currents_list_accumulated = zeros(Float64, steps+1, length(bonds))
density_correlations_accumulated = zeros(Float64, N, N)


K_avg = zeros(Int, steps, 9)
n_avg = zeros(Float64, steps+1, N)
avg_currents = zeros(Float64, steps+1, length(bonds))
avg_dd_correlations = zeros(Float64, N, N)
t_list = nothing
parameters = nothing

num_processes = 10
for run_idx in 1:num_processes

    filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))"

    data = load_from_hdf5(filename * "_run$(run_idx)" * ".h5")

    K_avg += data[:K_avg]
    n_avg += data[:n_avg]
    avg_currents += data[:avg_currents]
    avg_dd_correlations += data[:avg_dd_correlations]

    if run_idx == 1
        t_list = data[:t_list]
        parameters = data[:params]
    end

end

final_data = Dict(
        :K_avg => K_avg ./ num_processes,
        :n_avg => n_avg ./ num_processes,
        :avg_currents => avg_currents ./ num_processes,
        :avg_dd_correlations => avg_dd_correlations ./ num_processes,
        :t_list => t_list,
        :params => parameters
    )

filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_processes*num_iterations)_$(string(drive_type))_$(string(initial_state))"
save_to_hdf5(final_data, filename * ".h5")

# delete files after combining
for run_idx in 1:num_processes
    filename = "data/ff_$(Nx)x$(Ny)_dt$(dt)_p$(p)_b$(b)_t0.0_steps$(steps)_trajectories$(num_iterations)_$(string(drive_type))_$(string(initial_state))_run$(run_idx).h5"
    rm(filename)
end
