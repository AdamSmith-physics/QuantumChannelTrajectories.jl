

function trajectory(hamiltonian::SparseMatrixCSC, ψ::Vector, num_iterations::Int, params::SimulationParameters)

    steps = params.steps
    Nx = params.Nx
    Ny = params.Ny
    p = params.p
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

    for run in 1:num_iterations

        K_list = zeros(Int, steps, 9)
        n_list = zeros(Float64, steps+1, N)
        currents_list = zeros(Float64, steps+1, length(bonds))

        n_list[1, :] = [n_expectation(ψ, n, N) for n in 1:N]
        currents_list[1, :] = [current_expectation(ψ, B, bond, N) for bond in bonds]

        for step in 1:steps

            ψ, info = exponentiate(hamiltonian, -im*dt, ψ; tol=1e-15, ishermitian=true, eager=true) 

            # Kraus operator for inflow
            K_in = pick_kraus(p, ψ, site_in, N)
            apply_kraus!(ψ, K_in, site_in, N; type=:inflow, drive_type=drive_type)
        


       
    end


    return nothing

end