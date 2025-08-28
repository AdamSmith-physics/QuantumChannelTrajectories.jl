export create_circuit
export create_circuit_Anna
export create_circuit_Adam
export circuit_order



function create_circuit(Nx::Int, Ny::Int, Order; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)
    
    N::Int = Nx * Ny

    n_lyrs = length(Order)
    circuit_layers = [spzeros(Complex{Float64}, 2^N, 2^N) for _ in 1:n_lyrs]
    
    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    
    for lyr in 1:n_lyrs
        n_indxs = length(Order[lyr])
        for indxs in 1:n_indxs
            x1, x2 = Order[lyr][indxs]
            nx::Int = rem(x1,Nx)!==0 ? rem(x1,Nx) : Nx
            ny::Int = rem(x1,Nx)!==0 ? floor(x1/Nx)+1 : floor(x1/Nx)
            
            if abs(x1-x2)==1
                local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)
                row_operator = kron( sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)) )
                circuit_layers[lyr] += kron( sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1))) )        
            end    
            
            if abs(x1-x2)==Nx
                local_operator = V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
                local_operator -= kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
                local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
                row_operator = kron( sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx), 2^(Nx-nx)) )
                circuit_layers[lyr] += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))) )            
            end
        end
    end
    row_operator = nothing
    local_operator = nothing

    return circuit_layers
end




function create_circuit_Anna(Nx::Int, Ny::Int; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)

    N::Int = Nx * Ny

    layer_1 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_2 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_3 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_4 = spzeros(Complex{Float64}, 2^N, 2^N)

    row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
    for ny in 1:Ny
        local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)
        
        ########### Constructing first layer of horizontal operators ###########
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
        start = isodd(ny) ? 1 : 2
        for nx in start : 2 : Nx - 1
            row_operator += kron( sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)) )
        end
        layer_1 += kron( sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1))) )        
                
        ########### Constructing second layer of horizontal operators ###########
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)    
        start = isodd(ny) ? 2 : 1
        for nx in start:2:Nx-1
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)),local_operator,sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)))
        end
        layer_2 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),row_operator,sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1))))
    end    

    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    local_operator = V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
    local_operator -= kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
    local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
    
    ########### Constructing first layer of vertical operators ###########
    for ny in 1:Ny-1   
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        start = isodd(ny) ? 1 : 2
        for nx in start:2:Nx
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx), 2^(Nx-nx)))
        end
        layer_3 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))))
    # end
    ########### Constructing second layer of vertical operators ###########
    # for ny in 2:2:Ny-1 
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        start = isodd(ny) ? 2 : 1
        for nx in start:2:Nx
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx), 2^(Nx-nx)))
        end
        layer_4 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))))            
    end
            
    row_operator = nothing
    local_operator = nothing
            
    return layer_1, layer_3, layer_2, layer_4
end



function create_circuit_Adam(Nx::Int, Ny::Int; B::Float64 = 0.0, V::Float64 = 0.0, fermions::Bool = false)

    N::Int = Nx * Ny

    layer_1 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_2 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_3 = spzeros(Complex{Float64}, 2^N, 2^N)
    layer_4 = spzeros(Complex{Float64}, 2^N, 2^N)

    row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
    for ny in 1:Ny
        local_operator = -exp(im*B*ny)*kron(Sigma_plus,Sigma_minus) - exp(-im*B*ny)*kron(Sigma_minus,Sigma_plus) + V*kron(density_operator,density_operator)

        # Constructing first layer of odd horizontal operators
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)
        for nx in 1:2:Nx-1
            row_operator += kron( sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)) )
        end
        layer_1 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1))) )        
                
        # Constructing third layer of even horizontal operators
        row_operator = spzeros(Complex{Float64}, 2^Nx, 2^Nx)    
        for nx in 2:2:Nx-1
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx-1), 2^(Nx-nx-1)))
        end
        layer_3 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - Nx - Nx*(ny-1)), 2^(N - Nx - Nx*(ny-1))) )
    end    

    fill_operator = fermions ? PauliZ : sparse(I, 2, 2)
    local_operator = V*kron(density_operator, sparse(I,2^(Nx-1),2^(Nx-1)), density_operator)
    local_operator -= kron(Sigma_plus, fill(fill_operator,Nx-1)... ,Sigma_minus) 
    local_operator -= kron(Sigma_minus, fill(fill_operator,Nx-1)..., Sigma_plus) 
    
    # creating second layer of odd vertical operators    
    for ny in 1:2:Ny-1   
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        for nx in 1:Nx
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx), 2^(Nx-nx)))
        end
        layer_2 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))), row_operator, sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))))
    end
    
    # creating forth layer of even vertical operators
    for ny in 2:2:Ny-1 
        row_operator = spzeros(Complex{Float64}, 2^(2*Nx), 2^(2*Nx))
        for nx in 1:Nx
            row_operator += kron(sparse(I, 2^(nx-1), 2^(nx-1)), local_operator, sparse(I, 2^(Nx-nx), 2^(Nx-nx)))
        end
        layer_4 += kron(sparse(I, 2^(Nx*(ny-1)), 2^(Nx*(ny-1))),row_operator, sparse(I, 2^(N - 2*Nx - Nx*(ny-1)), 2^(N - 2*Nx - Nx*(ny-1))))            
    end
            
    row_operator = nothing
    local_operator = nothing
            
    return layer_1, layer_2, layer_3, layer_4
end


function circuit_order(Nx::Int,Ny::Int; type="Anna")
    
    layers=[[] for _ in 1:4]
    
    if type=="Anna"
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && (col+row)%2==0 ? push!(layers[1], (index,index+1)) : nothing    
                row<Ny && (col+row)%2==0 ? push!(layers[3], (index,index+Nx)) : nothing    
                col<Nx && (col+row)%2==1 ? push!(layers[2], (index,index+1)) : nothing    
                row<Ny && (col+row)%2==1 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    elseif type=="Adam"
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && col%2==1 ? push!(layers[1], (index,index+1)) : nothing    
                col<Nx && col%2==0 ? push!(layers[3], (index,index+1)) : nothing    
                row<Ny && row%2==1 ? push!(layers[2], (index,index+Nx)) : nothing    
                row<Ny && row%2==0 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    else   
        for row in 1:Ny
            for col in 1:Nx
                index = col+(row-1)*Nx
                col<Nx && col%2==1 ? push!(layers[1], (index,index+1)) : nothing    
                col<Nx && col%2==0 ? push!(layers[2], (index,index+1)) : nothing    
                row<Ny && row%2==1 ? push!(layers[3], (index,index+Nx)) : nothing    
                row<Ny && row%2==0 ? push!(layers[4], (index,index+Nx)) : nothing     
            end
        end
    end
    
    return layers
end
