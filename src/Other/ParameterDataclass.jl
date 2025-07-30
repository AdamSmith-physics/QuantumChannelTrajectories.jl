export SimulationParameters
export to_dict
export from_dict

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
end

# Constructor with default values
function SimulationParameters(; steps, Nx, Ny, p, bonds, site_in, site_out, dt, B,
                            drive_type = :current, 
                            initial_state = :random)
    return SimulationParameters(steps, Nx, Ny, dt, p, B, bonds, site_in, site_out, drive_type, initial_state)
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
        :initial_state => string(params.initial_state)
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
        initial_state = Symbol(dict[:initial_state])
    )
end