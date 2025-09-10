export SimulationParameters
export to_dict
export from_dict
export get_t_list

struct SimulationParameters
    steps::Int
    Nx::Int
    Ny::Int
    dt::Float64
    p::Float64
    B::Float64
    bonds::Vector{Tuple{Int,Int}}  # List of bonds as pairs of sites
    site_in::Int
    site_out::Int
    drive_type::Symbol  # :current, :dephasing
    initial_state::Symbol  # :checkerboard, :empty, :random, :custom
    even_parity::Bool  # Optional parameter for even parity
    pinned_corners::Bool  # Optional parameter for pinned corners
    single_shot::Bool  # Optional parameter for single shot measurements
    trotter_evolution::Bool  # Optional parameter for single shot measurements
    initialization::Any  # Optional parameter for initializing a density pattern
end

# Constructor with default values
function SimulationParameters(; steps, Nx, Ny, dt, p, B, bonds, site_in, site_out,
                            drive_type = :current, 
                            initial_state = :random,
                            even_parity = false,
                            pinned_corners = false,
                            single_shot = false,
                            trotter_evolution = false,
                            initialization = nothing,
                            )

    return SimulationParameters(steps, Nx, Ny, dt, p, B, bonds, site_in, site_out, drive_type, initial_state, even_parity, pinned_corners, single_shot, trotter_evolution, initialization)
end


function to_dict(params::SimulationParameters)
    return Dict(
        :steps => params.steps,
        :Nx => params.Nx,
        :Ny => params.Ny,
        :dt => params.dt,
        :p => params.p,
        :B => params.B,
        :bonds => params.bonds,
        :site_in => params.site_in,
        :site_out => params.site_out,
        :drive_type => string(params.drive_type),
        :initial_state => string(params.initial_state),
        :even_parity => params.even_parity,
        :pinned_corners => params.pinned_corners,
        :single_shot => params.single_shot,
        :trotter_evolution => params.trotter_evolution,
        :initialization => params.initialization
        )
end


function from_dict(dict::Dict)
    return SimulationParameters(
        steps = dict[:steps],
        Nx = dict[:Nx],
        Ny = dict[:Ny],
        dt = dict[:dt],
        p = dict[:p],
        B = dict[:B],
        bonds = dict[:bonds],
        site_in = dict[:site_in],
        site_out = dict[:site_out],
        drive_type = Symbol(dict[:drive_type]),
        initial_state = Symbol(dict[:initial_state]),
        even_parity = dict[:even_parity],
        pinned_corners = dict[:pinned_corners],
        single_shot = dict[:single_shot],
        trotter_evolution = dict[:trotter_evolution],
        initialization = dict[:initialization]
    )
end


function get_t_list(parameters::SimulationParameters)
    return [LinRange(0, parameters.steps * parameters.dt, parameters.steps + 1);]
end