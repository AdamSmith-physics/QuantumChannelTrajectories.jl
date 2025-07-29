using Revise
using BenchmarkTools
using LinearAlgebra

using QuantumChannelTrajectories


Nx = 4
Ny = 4

N = Nx * Ny

ψ = rand(Complex{Float64}, 2^(Nx * Ny))
normalize!(ψ)

kraus = pick_kraus(0.5,ψ, 1, Nx * Ny)

bonds = get_bonds(Nx, Ny, 1, 2)

x = SimulationParameters(
    steps=100, 
    Nx=Nx, 
    Ny=Ny, 
    dt=0.1,
    p=0.5, 
    B=1.0,
    bonds=get_bonds(Nx, Ny, 1, 2), 
    site_in=1, 
    site_out=2, 
    drive_type="current", 
    initial_state="random")

x.drive_type
`