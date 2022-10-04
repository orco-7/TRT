include("lib.jl")
using Plots

function ppf()
    #function new_lattice(x_min::Int, x_max::Int, y_min::Int, y_max::Int, resolution, Nt, u_l, L_l, ν_l, rho_l, nu)
    L = new_lattice(0, 5, 0, 52, 10, 100, 0.1, 1, 0.1, 1, 3/10)#10^(-5))
    L = initialize_lattice!(L)
    

    anim = @animate for i in 1:L.Nt
    
        L.ρ, L.u, L.v = @inbounds compute_macro(L.f, L.cx, L.cy)
        if i > 1
            L.u_right[:, 1] = L.u[:, end]
            #L.u_right[:, 2] = L.v[:, end]
            L.u_left[:, 1] = L.u[:, 1]
            #L.u_left[:, 2] = L.v[:, 1]
            L.ρ_left[:] = L.ρ[:, 1]
            L.ρ_right[:] = L.ρ[:, end]
        end

        L.v = zeros(size(L.u))
        L.fₑ = @inbounds equilibrium(L.u, L.v, L.cx, L.cy, L.w, L.ρ)
        #println(L.f)
        #readline()
        L.f⁺, L.f⁻, L.fₑ⁺, L.fₑ⁻ = @inbounds f_plus_minus(L.f, L.fₑ)
        #println(L.fₑ⁺)
        #readline()
        #println(L.f⁺)
        #readline()
    
        println("Iteration n°: ", i)
        L.f = @inbounds collision_bgk(L.f, L.ω⁺, L.fₑ)
        println("collision done")
        L.f = @inbounds streaming(L.f)
        println("streaming done")

        #L.f = bounce_back(L.f, L.boundary)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_s(L.u_bottom, L.u, L.v, L.ρ, L.f)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_n(L.u_bottom, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = @inbounds pressure_e(L.ρ_right, L.u_right, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = @inbounds pressure_w(L.ρ_left, L.u_left, L.u, L.v, L.ρ, L.f)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_nw(L.u, L.v, L.ρ, L.f)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_ne(L.u, L.v, L.ρ, L.f)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_sw(L.u, L.v, L.ρ, L.f)
        #L.f, L.ρ, L.u, L.v = @inbounds velocity_se(L.u, L.v, L.ρ, L.f)
        
        bounce_back(L.f, L.boundary)
        println("bounce back done")
        #L.f = @inbounds periodic_boundary(L.f)
        println("periodic boundary done")

        #### artfificially fixing it 
        #for i in 1:L.Ny, j in 1:L.Nx
        #    if L.u[i, j] > 10 || L.u[i, j] < 0  
        #        L.u[i, j] = 1.5
        #    end

        #    if L.ρ[i, j] > 10 || L.ρ[i, j] < 0
        #        L.ρ[i, j] = 1
        #    end
        #end

        println("plotting u...")
        p1 = heatmap(L.u, ylabel = "u")
        #println("plotting v...")
        #p2 = heatmap(L.v, ylabel = "v")
        println("plotting ρ...")
        p3 = heatmap(L.ρ, ylabel = "ρ")
        title_ = string("Time-step n°: ", i)
        plot(p1, p3, layout = (2, 1), title = title_)
        println("plotting done")

    end

    gif(anim, "test.gif", fps = 10)
    ############ OKAY
end


### this initializes the lattice with the proper boundary conditions for a Poiseulle Flow with 
### Re = 10
### note that for the velocities the first coloumn is u, the second one is v 
function initialize_lattice!(L::Lattice)
    
    L.boundary = add_north_wall!(L.boundary, L.Nx)
    L.boundary = add_south_wall!(L.boundary, L.Nx, L.Ny)
    L.ω⁺ = 0.6
    L.ρ_left = 1.05 .* ones(L.Ny)
    L.ρ_right = ones(L.Ny)
    L.u_left = 0.1*ones(L.Ny, 2)
    L.u_left[:, 2] = zeros(L.Ny)
    for i in 1:L.Ny
        L.u_left[i, 1] = 4L.uₗ * (i/L.Ny - (i/L.Ny)^2)
    end
    L.u_right = 0.1*ones(L.Ny, 2)
    L.u_right[:, 2] = zeros(L.Ny)
    L.u_top = zeros(L.Nx, 2)
    L.u_bottom = zeros(L.Nx, 2)

    L.fₑ = equilibrium(L.u, L.v, L.cx, L.cy, L.w, L.ρ)

    L.f⁺, L.f⁻, L.fₑ⁺, L.fₑ⁻ = f_plus_minus(L.f, L.fₑ)

    L.f = L.fₑ
    #println(L.f)
    #readline()
    return L
end
