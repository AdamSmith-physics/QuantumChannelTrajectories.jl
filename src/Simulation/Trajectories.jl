export trajectory
export run_trajectories


function trajectory(hamiltonian::SparseMatrixCSC, ψ_init::Vector, parameters::SimulationParameters)
    
    steps = parameters.steps
    Nx = parameters.Nx
    Ny = parameters.Ny
    dt = parameters.dt
    p = parameters.p
    B = parameters.B
    bonds = parameters.bonds
    site_in = parameters.site_in
    site_out = parameters.site_out
    drive_type = parameters.drive_type
    initial_state = parameters.initial_state

    N = Nx * Ny

    if initial_state == :random
        ψ = random_state(Nx, Ny)
    else
        ψ = copy(ψ_init)
    end

    K_list = zeros(Int, steps, 9)
    n_list = zeros(Float64, steps+1, N)
    currents_list = zeros(Float64, steps+1, length(bonds))
    density_correlations = zeros(Float64, N, N)

    n_list[1, :] = [n_expectation(ψ, n, N) for n in 1:N]
    currents_list[1, :] = [current_expectation(ψ, B, bond, Nx, Ny) for bond in bonds]

    for step in 1:steps

        ψ, info = exponentiate(hamiltonian, -im*dt, ψ; tol=1e-14, ishermitian=true, eager=true) 

        # Kraus operator for inflow
        K_in = pick_kraus(p, ψ, site_in, N)
        apply_kraus!(ψ, K_in, site_in, N; type=:inflow, drive_type=drive_type)

        K_out = pick_kraus(p, ψ, site_out, N)
        apply_kraus!(ψ, K_out, site_out, N; type=:outflow, drive_type=drive_type)

        K_list[step, K_in + 3*K_out + 1] = 1
        n_list[step+1, :] = [n_expectation(ψ, n, N) for n in 1:N]
        currents_list[step+1, :] = [current_expectation(ψ, B, bond, Nx, Ny) for bond in bonds]
    
        if step == steps[end]
            # Add density correlations!
        end

    end

    return K_list, n_list, currents_list, density_correlations

end


function run_trajectories(hamiltonian::SparseMatrixCSC, ψ_init::Vector, num_iterations::Int, parameters::SimulationParameters)

    steps = parameters.steps
    Nx = parameters.Nx
    Ny = parameters.Ny
    bonds = parameters.bonds

    N = Nx * Ny

    K_list_accumulated = zeros(Int, steps, 9)
    n_list_accumulated = zeros(Float64, steps+1, N)
    currents_list_accumulated = zeros(Float64, steps+1, length(bonds))
    density_correlations_accumulated = zeros(Float64, N, N)

    prog = Progress(num_iterations; desc="Running trajectories...", showspeed=true)
    for run in 1:num_iterations

        K_list, n_list, currents_list, density_correlations = trajectory(hamiltonian, ψ_init, parameters)

        K_list_accumulated += K_list
        n_list_accumulated += n_list
        currents_list_accumulated += currents_list
        density_correlations_accumulated += density_correlations

        next!(prog; showvalues = [("Completed trajectories", "$run / $num_iterations")])

    end

    data = Dict(
        :K_list => K_list_accumulated,
        :n_list => n_list_accumulated,
        :currents_list => currents_list_accumulated,
        :density_correlations => density_correlations_accumulated,
        :params => to_dict(parameters)
    )

    return data

end