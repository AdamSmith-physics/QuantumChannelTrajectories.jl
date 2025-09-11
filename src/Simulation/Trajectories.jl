export trajectory
export run_trajectories


function trajectory(hamiltonian, ψ_init::Vector, fermions::Bool, parameters::SimulationParameters)
    
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
    # optional parameters
    even_parity = parameters.even_parity
    pinned_corners = parameters.pinned_corners
    single_shot = parameters.single_shot
    trotter_evolution = parameters.trotter_evolution
    initialization = parameters.initialization

    N = Nx * Ny

    if initial_state == :random
        ψ = random_state(Nx, Ny; even_parity = even_parity, pinned_corners = pinned_corners, site_in = site_in, site_out = site_out, initialization=initialization)
    else
        ψ = copy(ψ_init)
    end

    ϕ = copy(ψ)
    speed_list = zeros(Float64, steps+1)

    K_list = zeros(Int, steps, 9)
    n_list = zeros(Float64, steps+1, N)
    currents_list = zeros(Float64, steps+1, length(bonds))
    dd_correlations = zeros(Float64, N, N)

    n_list[1,:] = [n_expectation(ψ, n, N; single_shot=single_shot) for n in 1:N]
    currents_list[1,:] = [current_expectation(ψ, B, bond, Nx, Ny, fermions; single_shot=single_shot) for bond in bonds]

    speed_list[1] = sqrt(1.0-abs(ϕ' * ψ)^2)/1.0
    # speed_list[1,2] = real( sqrt( ψ' * hamiltonian * hamiltonian * ψ - ψ' * hamiltonian * ψ))
    # prog = Progress(steps; dt=0.1, desc="Running trajectory...", showspeed=true)
    for step in 1:steps

        ## Evolve system by full or trotterized hamiltonian
        if trotter_evolution
            for layer in 1:4
                ψ, info = exponentiate(hamiltonian[layer], -im*dt, ψ; tol=1e-14, ishermitian=true, eager=true) 
            end
        else
            ψ, info = exponentiate(hamiltonian, -im*dt, ψ; tol=1e-14, ishermitian=true, eager=true) 
        end

        # Kraus operator for inflow
        K_in = pick_kraus(p, ψ, site_in, N)
        ψ = apply_kraus(ψ, K_in, site_in, N; type=:inflow, drive_type=drive_type, fermions)

        K_out = pick_kraus(p, ψ, site_out, N)
        ψ = apply_kraus(ψ, K_out, site_out, N; type=:outflow, drive_type=drive_type, fermions)

        K_list[step, K_in + 3*K_out + 1] = 1
        n_list[step+1, :] = [n_expectation(ψ, n, N; single_shot=single_shot) for n in 1:N]
        currents_list[step+1, :] = [current_expectation(ψ, B, bond, Nx, Ny, fermions; single_shot=single_shot) for bond in bonds]
        
        pure_dist = sqrt(1.0-abs(ϕ' * ψ)^2)/1.0
        speed_list[step+1] = pure_dist        
        # energy = ψ' * hamiltonian * ψ
        # enegry_2 = ψ' * hamiltonian * hamiltonian * ψ
        # energy_var = real( sqrt(enegry_2 - energy^2))
        # speed_list[step+1, 2] = energy_var
        
        if step == steps[end]
            dd_correlations = density_correlations(ψ, Nx, Ny)
        end

        ϕ = copy(ψ)

    end

    return K_list, n_list, currents_list, dd_correlations, speed_list

end


function run_trajectories(hamiltonian, ψ_init::Vector, num_iterations::Int, fermions::Bool, parameters::SimulationParameters; eager_saving::Bool = false, filename::String = "")

    if eager_saving && filename == ""
        error("Filename must be provided for eager saving.")
    end

    steps = parameters.steps
    Nx = parameters.Nx
    Ny = parameters.Ny
    bonds = parameters.bonds

    N = Nx * Ny

    K_list_accumulated = zeros(Int, steps, 9)
    n_list_accumulated = zeros(Float64, steps+1, N)
    n_sq_list_accumulated = zeros(Float64, steps+1, N)
    currents_list_accumulated = zeros(Float64, steps+1, length(bonds))
    currents_sq_list_accumulated = zeros(Float64, steps+1, length(bonds))
    dd_correlations_accumulated = zeros(Float64, N, N)
    
    speed_list_accumulated = zeros(Float64, steps+1)

    for run in 1:num_iterations
        
        t_start = time()
        
        println("Running trajectory $run / $num_iterations")
        flush(stdout)
        
        # K_list, n_list, currents_list, dd_correlations = trajectory(hamiltonian, ψ_init, fermions, parameters)
        K_list, n_list, currents_list, dd_correlations, speed_list = trajectory(hamiltonian, ψ_init, fermions, parameters)
        
        K_list_accumulated .+= K_list
        n_list_accumulated .+= n_list
        n_sq_list_accumulated .+= n_list.^2
        currents_list_accumulated .+= currents_list
        currents_sq_list_accumulated .+= currents_list.^2
        dd_correlations_accumulated .+= dd_correlations
        
        speed_list_accumulated .+= speed_list

        if eager_saving
            params_dict = to_dict(parameters)
            filter!(kv -> kv[2] !== nothing, params_dict)    
            
            data = Dict(
                :K_avg => K_list_accumulated,
                :n_avg => n_list_accumulated,
                :n_sq_avg => n_sq_list_accumulated,
                :avg_currents => currents_list_accumulated,
                :currents_sq_avg => currents_sq_list_accumulated,
                :avg_dd_correlations => dd_correlations_accumulated,
                :t_list => get_t_list(parameters),
                :params => params_dict,         # :params => to_dict(parameters),
                :completed_trajectories => run,
                :speed_avg => speed_list_accumulated,
            )
            save_to_hdf5(data, filename)
        end

        t_end = time()
        
        _print_update(t_start, t_end, run, num_iterations)

    end

    if eager_saving
        return nothing
    end
    
    params_dict = to_dict(parameters)
    filter!(kv -> kv[2] !== nothing, params_dict)

    data = Dict(
        :K_avg => K_list_accumulated,
        :n_avg => n_list_accumulated,
        :n_sq_avg => n_sq_list_accumulated,
        :avg_currents => currents_list_accumulated,
        :currents_sq_avg => currents_sq_list_accumulated,
        :avg_dd_correlations => dd_correlations_accumulated,
        :t_list => get_t_list(parameters),
        :params => params_dict,         # :params => to_dict(parameters),
        :completed_trajectories => num_iterations,
        :speed_avg => speed_list_accumulated,

    )

    return data

end


function _print_update(t_start::Float64, t_end::Float64, run::Int, num_iterations::Int)
    elapsed_time = t_end - t_start

    println("Completed trajectory $run / $num_iterations")
    println("Elapsed time for trajectory $run: $(elapsed_time) seconds")
    
    # remaining time in hours:minutes:seconds (e.g. 01:05:23)
    remaining_time = (num_iterations - run) * elapsed_time
    hours = Int(div(remaining_time, 3600))
    minutes = Int(div(remaining_time % 3600, 60))
    seconds = Int(round(remaining_time % 60))
    # remaining time in hours:minutes:seconds (e.g. 01:05:23)
    println("Estimated remaining time: $(lpad(hours, 2, '0')):$(lpad(minutes, 2, '0')):$(lpad(seconds, 2, '0'))")
    
    println("")
    flush(stdout)
end