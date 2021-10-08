function oneD_Propagation_of_acoustic_waves_of_air()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    Δt::Float64 = 5.e-6
    pseudoMaxIter::UInt32 = 5

    time_end::Float64 = 2.3e-3
    save_time::Float64 = 2.3e-3
    save_iteration::UInt32 = 10000

    discretizationScheme = "central"
    gravity = [0.0 0.0 0.0]

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    initial_p(x::Float64,y::Float64) = 1.e5
    initial_u(x::Float64,y::Float64) = 1.0
    initial_v(x::Float64,y::Float64) = 0.0
    initial_T(x::Float64,y::Float64) = 300.0
    initial_Y(x::Float64,y::Float64) = 0.0

    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "function"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    function left_U_BCValue(t::Float64)
        u₀ = 1.0
        Δu₀ = 0.01*u₀
        f = 2000.0
        return u₀ + Δu₀ * sin(2.0*pi*f*t)
    end

    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, Δt, pseudoMaxIter,
    save_time, save_iteration,time_end,
    discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end


function oneD_Propagation_of_acoustic_waves_of_water()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    Δt::Float64 = 1.e-6
    pseudoMaxIter::UInt32 = 5

    time_end::Float64 = 6.e-4
    save_time::Float64 = 6.e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    discretizationScheme = "central"
    gravity = [0.0 0.0 0.0]

    push!(p∞, 621780000.0)
    push!(cᵥ, 3610.0)
    push!(γ, 1.19)
    push!(b, 6.7212e-4)
    push!(q, -1177788.0)

    push!(p∞, 621780000.0)
    push!(cᵥ, 3610.0)
    push!(γ, 1.19)
    push!(b, 6.7212e-4)
    push!(q, -1177788.0)

    initial_p(x::Float64,y::Float64) = 1.e5
    initial_u(x::Float64,y::Float64) = 1.0
    initial_v(x::Float64,y::Float64) = 0.0
    initial_T(x::Float64,y::Float64) = 300.0
    initial_Y(x::Float64,y::Float64) = 0.0

    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "function"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    function left_U_BCValue(t::Float64)
        u₀ = 1.0
        Δu₀ = 0.01*u₀
        f = 6000.0
        return u₀ + Δu₀ * sin(2.0*pi*f*t)
    end

    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, Δt, pseudoMaxIter,
    save_time, save_iteration,time_end,
    discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end








function low_mach_number_riemann_problem()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [5.e-4]
    Δt_iters::Vector{UInt32} = [1000000]
    pseudoMaxIter::Vector{UInt32} = [10]
    
    time_end::Float64 = 0.01
    save_time::Float64 = 0.01
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 0.5
            return 10000.0
        else
            return 10000.85
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.5
            return 0.2
        else
            return 0.202
        end
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.5
            return 1.3888888888888888888888888888889
        else
            return 1.3890069444444444444444444444444
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if x < 0.5
            return 0.0
        else
            return 0.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "zeroGradient"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end







function gas_gas_subsonic_shock_tube()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-5]
    Δt_iters::Vector{UInt32} = [1000000]
    pseudoMaxIter::Vector{UInt32} = [10]
    
    time_end::Float64 = 8.e-4
    save_time::Float64 = 8.e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 288.93)
    push!(γ, 1.66)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 0.5
            return 2.e5
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.5
            return 0.0
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.5
            return 293.78201579487867027490913594
        else
            return 289.35185185185185185185185185185
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if x < 0.5
            return 1.0
        else
            return 0.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "zeroGradient"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end





function gas_gas_transonic_shock_tube()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [5.e-6]
    Δt_iters::Vector{UInt32} = [1000000]
    pseudoMaxIter::Vector{UInt32} = [10]
    
    time_end::Float64 = 6.e-4
    save_time::Float64 = 6.e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 288.93)
    push!(γ, 1.66)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 0.5
            return 5.e5
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.5
            return 200.0
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.5
            return 293.94669181270091168201390563503
        else
            return 289.35185185185185185185185185185
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if x < 0.5
            return 1.0
        else
            return 0.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "zeroGradient"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end





function Shock_wave_interaction()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [5.e-5, 1.e-4]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [10, 10]
    
    time_end::Float64 = 0.038#0.016 #
    save_time::Float64 = 0.038
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 0.1
            return 1000.0
        elseif 0.1 < x < 0.9
            return 0.01
        else
            return 100.0
        end
    end
    function initial_u(x::Float64,y::Float64)
        return 0.0
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.1
            return 3.4722222222222222222222222222222
        elseif 0.1 < x < 0.9
            return 3.4722222222222222222222222222222e-5
        else
            return 0.34722222222222222222222222222222
        end
    end
    function initial_Y(x::Float64,y::Float64)
        return 0.0
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 0.0
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "fixedValue"
    right_U_BCValue = 0.0
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, nothing, nothing, nothing,
    nothing, right_U_BCValue, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end




function Lax_test_case()

    Nx::UInt32 = 50
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 20.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [5.e-5, 1.e-4]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 0.01
    save_time::Float64 = 0.01
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 10.0
            return 3.528e5
        else
            return 5.71e4
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 10.0
            return 200.727
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 10.0
            return 2752.8089887640449438202247191011
        else
            return 396.52777777777777777777777777778
        end
    end
    function initial_Y(x::Float64,y::Float64)
        return 0.0
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "zeroGradient"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end



function Mac_3_test_case()

    Nx::UInt32 = 50
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 20.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-5, 5.e-5]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 0.006
    save_time::Float64 = 0.006
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "1st"
    spatial_discretizationScheme = "upwind"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    
    function initial_p(x::Float64,y::Float64)
        if x < 10.0
            return 1.0333e6
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 10.0
            return 290.93
        else
            return 1122.61
        end
    end
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 10.0
            return 930.21706565263734047763085875609
        else
            return 347.22222222222222222222222222222
        end
    end
    function initial_Y(x::Float64,y::Float64)
        return 0.0
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "zeroGradient"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    right_P_BCtype = "zeroGradient"
    right_U_BCtype = "zeroGradient"
    right_V_BCtype = "zeroGradient"
    right_T_BCtype = "zeroGradient"
    right_Y_BCtype = "zeroGradient"

    bottom_P_BCtype = "zeroGradient"
    bottom_U_BCtype = "slip"
    bottom_V_BCtype = "slip"
    bottom_T_BCtype = "zeroGradient"
    bottom_Y_BCtype = "zeroGradient"

    top_P_BCtype = "zeroGradient"
    top_U_BCtype = "slip"
    top_V_BCtype = "slip"
    top_T_BCtype = "zeroGradient"
    top_Y_BCtype = "zeroGradient"
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end





