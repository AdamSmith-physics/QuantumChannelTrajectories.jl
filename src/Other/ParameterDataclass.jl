export SimulationParameters

struct SimulationParameters
    steps::Int
    Nx::Int
    Ny::Int
    dt::Float64
    p::Float64
    B::Float64
    bonds::Vector{Tuple{Int, Int}}
    site_in::Int
    site_out::Int
    drive_type::String  # "current", "dephasing"
    initial_state::String  # "checkerboard", "empty", "random", "custom"
end

# Constructor with default values
function SimulationParameters(; steps, Nx, Ny, p, bonds, site_in, site_out, dt, B,
                            drive_type::String = "current", 
                            initial_state::String = "random")
    return SimulationParameters(steps, Nx, Ny, dt, p, B, bonds, site_in, site_out, drive_type, initial_state)
end