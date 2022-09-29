include("lib.jl")
using Plots

function pcf()
    #function new_lattice(x_min::Int, x_max::Int, y_min::Int, y_max::Int, resolution, Nt, u_l, L_l, ν_l, rho_l, nu)
    L = new_lattice(0, 5, 0, 2, 100, 100, 1, 1, 0.1, 1, 3/10)
    L = initialize_lattice!(L)
    

    anim = @animate for i in 1:L.Nt
    
        L.ρ, L.u, L.v = compute_macro(L.f, L.cx, L.cy)
        L.fₑ = equilibrium(L.u, L.v, L.cx, L.cy, L.w, L.ρ)
        #println(L.f)
        #readline()
        L.f⁺, L.f⁻, L.fₑ⁺, L.fₑ⁻ = f_plus_minus(L.f, L.fₑ)
        #println(L.fₑ⁺)
        #readline()
        #println(L.f⁺)
        #readline()
    
        println("Iteration n°: ", i)
        L.f = collision(L.f, L.f⁺, L.fₑ⁺, L.ω⁺, L.ω⁻, L.f⁻, L.fₑ⁻ )
        println("collision done")
        L.f = streaming(L.f)
        println("streaming done")

        #L.f = bounce_back(L.f, L.boundary)
        L.f, L.ρ, L.u, L.v = velocity_s(L.u_bottom, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = velocity_n(L.u_bottom, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = pressure_e(L.ρ_right, L.u_right, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = pressure_w(L.ρ_left, L.u_left, L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = velocity_nw(L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = velocity_ne(L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = velocity_sw(L.u, L.v, L.ρ, L.f)
        L.f, L.ρ, L.u, L.v = velocity_se(L.u, L.v, L.ρ, L.f)
        println("bounce back done")

        L.f = periodic_boundary(L.f)
        println("periodic boundary done")

        #### artfificially fixing it 
        for i in 1:L.Ny, j in 1:L.Nx
            if L.u[i, j] > 10 || L.u[i, j] < 0  
                L.u[i, j] = 1.5
            end

            if L.ρ[i, j] > 4 || L.ρ[i, j] < 0
                L.ρ[i, j] = 1
            end
        end

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

    gif(anim, "test_pcf.gif", fps = 10)
end

function initialize_lattice!(L::Lattice)
    
    #L.boundary = add_north_wall!(L.boundary, L.Nx)
    #L.boundary = add_south_wall!(L.boundary, L.Nx, L.Ny)
    #println(L.Ny)
    #println(L.Nx)
    L.ρ_left =  ones(L.Ny)
    L.ρ_right = ones(L.Ny)
    #println("size ρ_left: ", size(L.ρ_left), ", size of ρ_right: ", size(L.ρ_right))
    L.u_left = ones(L.Ny, 2)
    L.u_right = ones(L.Ny, 2)
    L.u_top = ones(L.Nx, 2)
    L.u_top[:, 2] = zeros(L.Nx)
    L.u_bottom = -1 .* ones(L.Nx, 2)
    L.u_bottom[:, 2] = zeros(L.Nx)

    L.fₑ = equilibrium(L.u, L.v, L.cx, L.cy, L.w, L.ρ)

    L.f⁺, L.f⁻, L.fₑ⁺, L.fₑ⁻ = f_plus_minus(L.f, L.fₑ)

    L.f = L.fₑ
    #println(L.f)
    #readline()
    return L
end
