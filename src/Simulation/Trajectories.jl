export trajectory

function trajectory(hamiltonian::SparseMatrixCSC, ψ::Vector, num_iterations::Int, params::SimulationParameters)

    BLAS.set_num_threads(1) 

    steps = params.steps
    Nx = params.Nx
    Ny = params.Ny
    dt = params.dt
    p = params.p
    B = params.B
    bonds = params.bonds
    site_in = params.site_in
    site_out = params.site_out
    drive_type = params.drive_type
    initial_state = params.initial_state

    N = Nx * Ny

    K_list_accumulated = zeros(Int, steps, 9)
    n_list_accumulated = zeros(Float64, steps+1, N)
    currents_list_accumulated = zeros(Float64, steps+1, length(bonds))
    density_correlations_accumulated = zeros(Float64, N, N)

    prog = Progress(num_iterations; desc="Running trajectories...", showspeed=true)
    for run in 1:num_iterations

        K_list = zeros(Int, steps, 9)
        n_list = zeros(Float64, steps+1, N)
        currents_list = zeros(Float64, steps+1, length(bonds))
        density_correlations = zeros(Float64, N, N)

        n_list[1, :] = [n_expectation(ψ, n, N) for n in 1:N]
        currents_list[1, :] = [current_expectation(ψ, B, bond, Nx, Ny) for bond in bonds]

        for step in 1:steps

            ψ, info = exponentiate(hamiltonian, -im*dt, ψ; tol=1e-15, ishermitian=true, eager=true) 

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
        :density_correlations => density_correlations_accumulated
    )

    return data

end