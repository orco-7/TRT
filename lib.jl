include("lattice.jl")


### in this file are contained all the core functions of the algorithm

##########################################################################################
##########################################################################################

### computes the macroscopic quantities and overwrites them in the 
### respective field of L 
function compute_macro(f, cx, cy)
    ρ = zeros(size(f)[1:2])
    u_ = zeros(size(f)[1:2])
    v_ = zeros(size(f)[1:2])
    temp_ρ = sum(f, dims = 3)
    ρ[:, :] = temp_ρ[:, :, 1]
    for i in 1:size(f)[3]
        u_ += cx[i] .* f[:, :, i] 
        v_ += cy[i] .* f[:, :, i] 
    end
    u = u_ ./ ρ
    v = v_ ./ ρ
    return ρ, u, v
    ######## OKAY
end

### computes the symmetric and antisymmetric parts of the distributions
function f_plus_minus(f, fₑ)
    ### resting population have 0 antisymmetric population 
    ### the other populations are computed in the for cycle  
    f⁺ = zeros(size(f))
    f⁻ = zeros(size(f))
    fₑ⁺ = zeros(size(f))
    fₑ⁻ = zeros(size(f))

    f⁺[:, :, 1] = f[:, :, 1]
    fₑ⁺[:, :, 1] = fₑ[:, :, 1]
    f⁻[:, :, 1] .= 0
    fₑ⁻[:, :, 1] .= 0
    ### in position i there's the index of its opposite velocity 
    ind = [1, 3, 2, 5, 4, 7, 6, 9, 8]

    for i in 2:size(f)[3]
        j = ind[i]
        f⁺[:, :, i] = (f[:, :, i] + f[:, :, j]) ./ 2
        f⁻[:, :, i] = (f[:, :, i] - f[:, :, j]) ./ 2
        fₑ⁺[:, :, i] = (fₑ[:, :, i] + fₑ[:, :, j]) ./ 2
        fₑ⁻[:, :, i] = (fₑ[:, :, i] - fₑ[:, :, j]) ./ 2
    end
    return f⁺, f⁻, fₑ⁺, fₑ⁻
    ########### OKAY
end


### computes the equilibrium distribution and assigns it to the proper 
### field in the Lattice struct 
function equilibrium(u, v, cx, cy, w, ρ)
    ### this two temp variables are used exclusively to not repeat the same operations multiple times 
    fₑ = zeros(size(u)[1], size(v)[2], 9)
    temp = 1.5.*(u.*u + v.*v)
    for i in 1:9
        temp2 = 3 .* (u.*cx[i] + v.*cy[i])
        fₑ[:, :, i] = w[i] .* ρ .* (1 .+ temp2 + 0.5.* temp2.^2 - temp)
    end
    return fₑ
    ######### OKAY 
end

#### this function computes collision with the BGK operator, it uses the same lattice with tau = ω⁺
function collision_bgk(f, ω, fₑ)
    f_temp = zeros(size(f))
    ω_ = 1 - ω
    for i in 1:size(f_temp)[3]
        f_temp[:, :, i] = f[:, :, i] .* ω_ .+ fₑ[:, :, i].*ω
    end
    return f_temp
end

### computes collision with the TRT operator 
function collision(f, f⁺, fₑ⁺, ω⁺, ω⁻, f⁻, fₑ⁻ )
    
    ### buffer to prevent memory corruption, shouldn't result in too much overhead 
    ### after benchmark I might change it 
    f_temp = f - ω⁺ .* (f⁺ - fₑ⁺) - ω⁻ .* (f⁻ - fₑ⁻)
    return f_temp
    ####### OKAY 
end

### computes streaming 
function streaming(f)
    ### buffer, same as before 
    temp_f = zeros(size(f))
    temp_f[:, :, 1] = f[:, :, 1]
    temp_f[:, 2:end, 2] = f[:, 1:end-1, 2] 
    temp_f[:, 1:end-1, 3] = f[:, 2:end, 3]
    temp_f[1:end-1, :, 4] = f[2:end, :, 4]
    temp_f[2:end, :, 5] = f[1:end-1, :, 5]
    temp_f[1:end-1, 2:end, 6] = f[2:end, 1:end-1, 6]
    temp_f[2:end, 1:end-1, 7] = f[1:end-1, 2:end, 7]
    temp_f[1:end-1, 1:end-1, 8] = f[2:end, 2:end, 8]
    temp_f[2:end, 2:end, 9] = f[1:end-1, 1:end-1, 9]
    return temp_f
    ######## OKAY
end

#############################################################################################################################
#############################################################################################################################

### Here are defined the functions for the Zou-He and bounce-back boundary conditions, this are essentialy Dirichlet boundary 
### conditions (no-slip) plus pressure and velocity boundary conditions at inlets and outlets (note that this is an 
### incompressible simulation so pressure is directly proportional to density and the density values will be used. 
### Pressure can be obtained from density via a factor of 1/cₛ²). Note that Zou-He boundary conditions should only be applied 
### only if at the boundary the velcity is non-zero, for stationary walls simple bounce back is more stable and more accurate 

### the location of the boundaries will be denoted with cardinal directions so for a velocity boundary condition
### on an upper boundary velocity_n() will be used, upper-left will use velovity_ne() and so on


### This is a first implementation and only supports walls at the extremities of the domain
### they will need to be changed to accept coordinates of a body as parameters, to implement once the algorithm actually works 

### Standard link-wise bounce-back
function bounce_back(f, boundary)
    ### in position i there's the index of its opposite velocity 
    ind::Array{Int, 1} = [1, 3, 2, 5, 4, 7, 6, 9, 8]
    f_temp = f

    for n in 2:size(boundary)[2]
        i::Int = boundary[1, n]
        j::Int = boundary[2, n]
        q::Int = boundary[3, n]

        q_ = ind[q]

        f[i, j, q] = f_temp[i, j, q_]
    end
    return f 
    ####### OKAY
end


### periodic boundary conditions 
function periodic_boundary(f)
    f[:, 1, 6] = f[:, end, 6]
    f[:, 1, 2] = f[:, end, 2] 
    f[:, 1, 9] = f[:, end, 9]

    f[:, end, 7] = f[:, 1, 7]
    f[:, end, 3] = f[:, 1, 3]
    f[:, end, 8] = f[:, 1, 8]
    return f  
end

### only this function will be commented, the others are basically the same 
function velocity_n(u_top,u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2
    ### the macroscopic velocities at the boundary are directly imposed 
    u[1, :] = u_top[:, 1]
    v[1, :] = u_top[:, 2]
    ### density is obtained form the known distributions, unknown distributions are a combination of velocity, density 
    ### and known distributions, see "On Pressure and Velocity Boundary Conditions ..." pg.4 for analitical proof 

    ### remainder of the velocity's notation 
    ###        8 --- 4 --- 6                                                 ###
    ###        |  \  |  /  |                                                 ###
    ###        3 --- 1 --- 2                                                 ###
    ###        |  /  |  \  |                                                 ###
    ###        7 --- 5 --- 9                                                 ###
    ρ[1, 2:end-1] = (f[1, 2:end-1, 1] + f[1, 2:end-1, 2] + f[1, 2:end-1, 3] + 2 .*(f[1, 2:end-1, 4] + 
                f[1, 2:end-1, 8] + f[1, 2:end-1, 6])) ./ (1 .+ u[1, 2:end-1])

    ### need to bemchmark to find out if it's less expensive to store the values needed in separate vectors to use in computations 

    f[1, 2:end-1, 5] = f[1, 2:end-1, 4] - k1 .* (ρ[1, 2:end-1] .* u[1, 2:end-1])
    f[1, 2:end-1, 7] = f[1, 2:end-1, 6] + k3 .* (f[1, 2:end-1, 2] - f[1, 2:end-1, 3]) + k2 .* (v[1, 2:end-1] .* ρ[1, 2:end-1]) - 
                            k3 .* (u[1, 2:end-1].*ρ[1, 2:end-1])
    f[1, 2:end-1, 9] = f[1, 2:end-1, 8] - k3 .* (f[1, 2:end-1, 2] - f[1, 2:end-1, 3]) - k2 .* (v[1, 2:end-1] .* ρ[1, 2:end-1]) + 
                            k3 .* (u[1, 2:end-1].*ρ[1, 2:end-1])
    
    
    return f, ρ, u, v

end

function velocity_s(u_bottom, u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2
    u[end, :] = u_bottom[:, 1]
    v[end, :] = u_bottom[:, 2]

    ρ[end, 2:end-1] = (f[end, 2:end-1, 1] + f[end, 2:end-1, 2] + f[end, 2:end-1, 3] + 2 .*(f[end, 2:end-1, 5] + 
                f[end, 2:end-1, 7] + f[end, 2:end-1, 9])) ./ (1 .+ u[end, 2:end-1])

    f[end, 2:end-1, 4] = f[end, 2:end-1, 5] + k1 .* (ρ[end, 2:end-1] .* u[end, 2:end-1])
    f[end, 2:end-1, 6] = f[end, 2:end-1, 7] - k3 .* (f[end, 2:end-1, 2] - f[end, 2:end-1, 3]) + k2 .*(v[end, 2:end-1] .* ρ[end, 2:end-1]) +
                                k3 .* (u[end, 2:end-1] .* ρ[end, 2:end-1])
    f[end, 2:end-1, 8] = f[end, 2:end-1, 9] + k3 .* (f[end, 2:end-1, 2] - f[end, 2:end-1, 3]) + k2 .*(v[end, 2:end-1] .* ρ[end, 2:end-1]) -
                                k3 .* (u[end, 2:end-1] .* ρ[end, 2:end-1])
    return f, ρ, u, v
end

function velocity_e(u_right, u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2
    u[:, end] = u_right[:, 1]
    v[:, end] = u_right[:, 2]

    ρ[2:end-1, end] = (f[2:end-1, end, 1] + f[2:end-1, end, 4] + f[2:end-1, end, 5] + 2 .* (f[2:end-1, end, 3] + f[2:end-1, end, 7] + 
                       f[2:end-1, end, 8])) / (1 .+ v[2:end-1, end]) 
    
    f[2:end-1, end, 3] = f[2:end-1, end, 2] - k1 .* (ρ[2:end-1, end] .* u[2:end-1, end])
    f[2:end-1, end, 7] = f[2:end-1, end, 6] - k3 .* (f[2:end-1, end, 4] + f[2:end-1, end, 5]) + k2 .* (u[2:end-1, end] .* ρ[2:end-1, end]) + 
                         k3 .* (v[2:end-1, end] .* ρ[2:end-1, end])
    f[2:end-1, end, 8] = f[2:end-1, end, 9] + k3 .* (f[2:end-1, end, 4] + L.f[2:end-1, end, 5]) + k2 .* (L.u[2:end-1, end] .* L.ρ[2:end-1, end]) - 
                         k3 .* (L.v[2:end-1, end] .* L.ρ[2:end-1, end])
    return f, ρ, u, v
end

function velocity_w(u_left, u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2
    u[:, 1] = u_left[:, 1]
    v[:, 1] = u_left[:, 2]

    ρ[2:end-1, 1] = (f[2:end-1, 1, 1] + f[2:end-1, 1, 4] + f[2:end-1, 1, 5] + 2 .* (f[2:end-1, 1, 2] + f[2:end-1, 1, 6] + 
                       f[2:end-1, 1, 9])) / (1 .+ v[2:end-1, 1]) 
    
    f[2:end-1, 1, 2] = f[2:end-1, 1, 3] + k1 .* (ρ[2:end-1, 1] .* u[2:end-1, 1])
    f[2:end-1, 1, 6] = f[2:end-1, 1, 7] + k3 .* (f[2:end-1, 1, 4] + f[2:end-1, 1, 5]) + k2 .* (u[2:end-1, 1] .* ρ[2:end-1, 1]) -
                             k3 .* (v[2:end-1, 1] .* ρ[2:end-1, 1])
    f[2:end-1, 1, 9] = f[2:end-1, 1, 8] - k3 .* (f[2:end-1, 1, 4] + f[2:end-1, 1, 5]) + k2 .* (u[2:end-1, 1] .* ρ[2:end-1, 1]) + 
                         k3 .* (v[2:end-1, 1] .* ρ[2:end-1, 1])
    return f, ρ, u, v
end

function pressure_w(ρ_left, u_left, u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2
    ### for pressure boundary conditions density along the boundary and the component of velocity tangential to it 
    ### are imposed 
    ρ[:, 1] = ρ_left[:]
    v[:, 1] = u_left[:, 2]

    ### now we can determine the normal component of the velocity similarly to how density is determined 

    u[2:end-1, 1] = (f[2:end-1, 1, 1] + f[2:end-1, 1, 4] + f[2:end-1, 1, 5] + 2 .* (f[2:end-1, 1, 2] + f[2:end-1, 1, 6] + 
                       f[2:end-1, 1, 9])) ./ ( 1 .+ v[2:end-1, 1])
    
    ### now the unknown distributions are computed in the same way as for the velocity bc
    f[2:end-1, 1, 2] = f[2:end-1, 1, 3] - k1 .* (ρ[2:end-1, 1] .* u[2:end-1, 1])
    f[2:end-1, 1, 6] = f[2:end-1, 1, 7] + k3 .* (f[2:end-1, 1, 4] + f[2:end-1, 1, 5]) + k2 .* (u[2:end-1, 1] .* ρ[2:end-1, 1]) -
                         k3 .* (v[2:end-1, 1] .* ρ[2:end-1, 1])
    f[2:end-1, 1, 9] = f[2:end-1, 1, 8] - k3 .* (f[2:end-1, 1, 4] + f[2:end-1, 1, 5]) + k2 .* (u[2:end-1, 1] .* ρ[2:end-1, 1]) + 
                         k3 .* (v[2:end-1, 1] .* ρ[2:end-1, 1])
    return f, ρ, u, v
end
### same function but for the left boundary 
function pressure_e(ρ_right, u_right, u, v, ρ, f)
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 2/3
    k2 = 1/6
    k3 = 1/2

    ρ[:, end] = ρ_right[:, 1]
    v[:, end] = u_right[:, 2]

    u[2:end-1, end] = (f[2:end-1, end, 1] + f[2:end-1, end, 4] + f[2:end-1, end, 5] + 2 .* (f[2:end-1, end, 3] + f[2:end-1, end, 7] + 
                       f[2:end-1, end, 8])) ./ (1 .+ v[2:end-1, end]) 
    
    f[2:end-1, end, 3] = f[2:end-1, end, 2] - k1 .* (ρ[2:end-1, end] .* u[2:end-1, end])
    f[2:end-1, end, 7] = f[2:end-1, end, 6] - k3 .* (f[2:end-1, end, 4] + f[2:end-1, end, 5]) + k2 .* (u[2:end-1, end] .* ρ[2:end-1, end]) + 
                         k3 .* (v[2:end-1, end] .* ρ[2:end-1, end])
    f[2:end-1, end, 8] = f[2:end-1, end, 9] + k3 .* (f[2:end-1, end, 4] + f[2:end-1, end, 5]) + k2 .* (u[2:end-1, end] .* ρ[2:end-1, end]) - 
                         k3 .* (v[2:end-1, end] .* ρ[2:end-1, end])
    return f, ρ, u, v
end


### now are implemented the boundary conditions for the corners 
function velocity_nw(u, v, ρ, f)
    ### velocity and density in the corner are simply interpoleted from the sourronding ones 
    u[1, 1] = u[2, 2]
    v[1, 1] = v[2, 2]
    ρ[1, 1] = ρ[2, 2]
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 0.5
    k2 = 1/6
    k3 = 2/3
    ### we now determine the populations "streaming" from the boundaries, for analytical proof refer to the same 
    ### paper as before and to "Consistent Boundary Conditions for 2D and 3D LBM ..; Chao-An Lin"
    f[1, 1, 7] = k1*(ρ[1, 1]*(1 - u[1, 1] + k2*v[1, 1]) - f[1, 1, 1] - 2*(f[1, 1, 3] + f[1, 1, 4] + f[1, 1, 8]))
    f[1, 1, 2] = f[1, 1, 3] + k3*ρ[1, 1]*u[1, 1]
    f[1, 1, 5] = f[1, 1, 4] - k3*ρ[1, 1]*v[1, 1]
    f[1, 1, 6] = f[1, 1, 7] + k2*ρ[1, 1]*(u[1, 1] + v[1, 1])
    f[1, 1, 9] = f[1, 1, 8] + k2*ρ[1, 1]*(u[1, 1] + v[1, 1])
    return f, ρ, u, v
end

function velocity_ne(u, v, ρ, f)
    u[1, end] = u[2, end-1]
    v[1, end] = v[2, end-1]
    ρ[1, end] = ρ[2, end-1]
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 0.5
    k2 = 1/6
    k3 = 2/3

    f[1, end, 9] = k1*(ρ[1, end] * (1 - u[1, end] + k2*v[1, end]) - f[1, end, 1] - 2*(f[1, end, 2] + f[1, end, 4] + f[1, end, 6]))
    f[1, end, 3] = f[1, end, 2] + k3*ρ[1, end]*u[1, end]
    f[1, end, 5] = f[1, end, 4] - k3*ρ[1, end]*v[1, end]
    f[1, end, 7] = f[1, end, 6] + k2*ρ[1, end]*(u[1, end] + v[1, end])
    f[1, end, 8] = f[1, end, 9] + k2*ρ[1, end]*(u[1, end] + v[1, end])
    return f, ρ, u, v
end

function velocity_sw(u, v, ρ, f)
    u[end, 1] = u[end-1, 2]
    v[end, 1] = u[end-1, 2]
    ρ[end, 1] = ρ[end-1, 2]
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 0.5
    k2 = 1/6
    k3 = 2/3

    f[end, 1, 8] = k1*(ρ[end, 1] * (1 - u[end, 1] + k2*v[end, 1]) - f[end, 1, 1] - 2*(f[end, 1, 3] + f[end, 1, 5] + f[end, 1, 7]))
    f[end, 1, 2] = f[end, 1, 3] + k3*ρ[end, 1]*u[end, 1]
    f[end, 1, 4] = f[end, 1, 5] - k3*ρ[end, 1]*v[end, 1]
    f[end, 1, 6] = f[end, 1, 7] + k2*ρ[end, 1]*(u[end, 1] + v[end, 1])
    f[end, 1, 9] = f[end, 1, 8] + k2*ρ[end, 1]*(u[end, 1] + v[end, 1])
    return f, ρ, u, v
end

function velocity_se(u, v, ρ, f)
    u[end, end] = u[end-1, end-1]
    v[end, end] = v[end-1, end-1]
    ρ[end, end] = ρ[end-1, end-1]
    ### this coefficients are used to slim down notation and to not repeat numerical division many times 
    k1 = 0.5
    k2 = 1/6
    k3 = 2/3

    f[end, end, 6] = k1*(ρ[end, end] * (1 - u[end, end] + k2*v[end, end]) - f[end, end, 1] - 2*(f[end, end, 2] + f[end, end, 5] + f[end, end, 9]))
    f[end, end, 4] = f[end, end, 5] - k3*ρ[end, end]*v[end, end]
    f[end, end, 3] = f[end, end, 2] + k3*ρ[end, end]*u[end, end]
    f[end, end, 7] = f[end, end, 6] + k2*ρ[end, 1]*(u[end, end] + v[end, end])
    f[end, end, 8] = f[end, end, 9] + k2*ρ[end, 1]*(u[end, end] + v[end, end])
    return f, ρ, u, v
end



### TODO: there are multiple inefficiencies in this version. A significant spped-up could be obtained by 
### using the @inbound macro, and by storing the needed values in the Zou-He functions in temporary vectors 
### in order to not read a coloumn/row from a matrix multiple times (this needs to be benchmarked against memory usage in
### a larger lattice). 

### All the variables will need to be statically typed to slim down compilation time and type inference


### There's also a significant bottolneck in plotting, need to check against memory allocations if storing the values of velocities
### and density in a (Ny, Nx, Nt) tensor is more convenient. 

### Also need to check a way to store permanently macroscopic variables in an efficient way. CSV doesn't seem to handle this well 
