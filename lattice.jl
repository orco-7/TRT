### this file defines the lattice struct and its building functions  
    ############################################################################
    ### this is a D2Q9 TRT implementation, the following scheme for numbering ##
    ### velocities will be used                                              ###
    ###        8 --- 4 --- 6                                                 ###
    ###        |  \  |  /  |                                                 ###
    ###        3 --- 1 --- 2                                                 ###
    ###        |  /  |  \  |                                                 ###
    ###        7 --- 5 --- 9                                                 ###
    ############################################################################
mutable struct Lattice
    #### grid parameters 
    Nx::Int
    Ny::Int
    Q::Int
    ### lattice parameters
    uₗ::Float64
    Lₗ::Float64
    νₗ::Float64
    Reₗ::Float64
    ρₗ::Float64
    cₛ::Float64
    ### number of iterations
    Nt::Int

    ### TRT parameters
    Λ::Float64
    ω⁺::Float64
    ω⁻::Float64

    #### velocities and relative weights in the same direction 
    cx::Array
    cy::Array
    w ::Array
    
    #### distributions
    f::Array{Float64, 3}
    f⁺::Array{Float64, 3}
    f⁻::Array{Float64, 3}
    fₑ::Array{Float64, 3}
    fₑ⁺::Array{Float64, 3}
    fₑ⁻::Array{Float64, 3}


    ### Boundary conditions
    ### this are (Ny/Nx, 2)ᵗ arrays, the first row contains 
    ### the x component of the velocity while the second one the 
    ### y component, ρ_right(left) is simply a (Ny, 1) Array  
    u_top::Array{Float64, 2}
    u_bottom::Array{Float64, 2}
    u_left::Array{Float64, 2}
    u_right::Array{Float64, 2}
    ρ_right::Array{Float64, 1}
    ρ_left::Array{Float64, 1}
    ### this initially has the same dimension as the matrix, when all the boundaries and bodies are defined 
    ### will store in each coloumn the the inidices of the nodes where bounce-back must be applied so it will be 
    ### an Array{Int, 2} 2xn_boundary 
    boundary

    ### macroscopic quantities 
    ρ
    u
    v
end


function new_lattice(x_min::Int, x_max::Int, y_min::Int, y_max::Int, resolution, Nt, u_l, L_l, ν_l, rho_l, nu)
    Nx = resolution*(x_max-x_min)
    Ny = resolution*(y_max-y_min) 

    Re_l = u_l*L_l/ν_l
    println("Reynolds number = ", Re_l)
    
    ### for stability reasons the "magic parameter" is fixed at 0.4, with stability analysis this could be changed 
    ### om_m is then linearly dipendent from lambda and omega plus 
    lam = 0.4
    ### ω⁺ is a function of viscosity of the phisical system, therefore 
    om_p = (3*nu + 0.5)^(-1)
    om_m = (0.5 + lam/(1/(om_p - 0.5)))
    cx = [0, 1, -1, 0, 0, 1, -1, -1, 1]
    cy = [0, 0, 0, 1, -1, 1, -1, 1, -1]
    w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]

    f_prot = zeros(Ny, Nx, 9)
    u_v = zeros(Ny, 2)
    u_h = zeros(Nx, 2)
    rho_prot = zeros(Ny)
    l_prot = zeros(Ny, Nx)

    L = Lattice(Nx, Ny, 9, u_l, L_l, ν_l, Re_l, rho_l, 1/sqrt(3), Nt, lam, om_p, om_m, cx, cy, w, f_prot, f_prot, f_prot, f_prot, f_prot, f_prot, u_h, u_h,
               u_v, u_v, rho_prot, rho_prot, zeros(3), ones(Ny, Nx), l_prot, l_prot)
    
    return L
end


##########################################################################################################################################
#################################### get new boundaries ################################################################################## 

function add_north_wall!(boundary, Nx)
    for i in 1:Nx, q in [4, 6, 8]
        boundary = hcat(boundary, [1, i, q])
    end
    return boundary 
end

function add_south_wall!(boundary, Nx, Ny)
    for i in 1:Nx, q in [5, 7, 9]
        boundary = hcat(boundary, [Ny, i, q])
    end
    return boundary 
end

function add_east_wall!(boundary, Nx, Ny)
    for i in 1:Ny, q in [2, 6, 9]
        boundary = hcat(boundary, [i, Nx, q])
    end
    return boundary 
end

function add_west_wall!(boundary, Ny)
    for i in 1:Ny, q in [3, 7, 8]
        boundary = hcat(boundary, [i, 1, q])
    end
    return boundary 
end