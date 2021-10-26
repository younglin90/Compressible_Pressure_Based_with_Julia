

function oneD_test1_transmission_at_fluid_interfaces()


    Nx::UInt32 = 500
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.5
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [2.e-6, 2.e-6]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 1.75e-3
    save_time::Float64 = 0.75e-3
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "linear"
    spatial_discretizationScheme_U = "vanleer"
    spatial_discretizationScheme_T = "minmod"
    spatial_discretizationScheme_Y = "minmod"

    gravity = [0.0 0.0 0.0]
    
    push!(p∞, 0.0)
    push!(cᵥ, 3115.6)
    push!(γ, 1.667)
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
    function initial_Y(x::Float64,y::Float64)
        if x < 1.0
            return 1.0
        else
            return 0.0
        end
    end

    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "function"
    left_V_BCtype = "zeroGradient"
    left_T_BCtype = "zeroGradient"
    left_Y_BCtype = "zeroGradient"
    
    function left_U_BCValue(t::Float64)
        u₀ = 1.0
        Δu₀ = 0.02*u₀
        f = 5000.0

        if t<1.0/f
            U = u₀ + Δu₀ * sin(pi*f*t)
            if U < u₀
                return u₀
            else
                return U
            end
        else
            return u₀ 
        end
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
    
    return Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
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




function oneD_test2_interface_advection_with_constant_velocity()

    Nx::UInt32 = 500
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [5.e-4, 5.e-4]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 0.7
    save_time::Float64 = 0.7
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "upwind"
    spatial_discretizationScheme_U = "upwind"
    spatial_discretizationScheme_T = "upwind"
    spatial_discretizationScheme_Y = "upwind"
    #spatial_discretizationScheme_Y = "minmod"
    #spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 3115.6)
    push!(γ, 1.667)
    push!(b, 0.0)
    push!(q, 0.0)


    initial_p(x::Float64,y::Float64) = 1.e5
    initial_u(x::Float64,y::Float64) = 1.0
    initial_v(x::Float64,y::Float64) = 0.0
    initial_T(x::Float64,y::Float64) = 300.0
    function initial_Y(x::Float64,y::Float64)
        if x < 0.1
            return 1.0
        else
            return 0.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 1.0
    left_V_BCtype = "fixedValue"
    left_V_BCValue = 0.0
    left_T_BCtype = "fixedValue"
    left_T_BCValue = 300.0
    left_Y_BCtype = "fixedValue"
    left_Y_BCValue = 1.0

    right_P_BCtype = "fixedValue"
    right_P_BCValue = 1.e5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, left_V_BCValue, left_T_BCValue, left_Y_BCValue,
    right_P_BCValue, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end







#=
function oneD_test3_air_helium_interface()

    Nx::UInt32 = 400
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 0.4
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-7, 1.e-7]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 2.e-4
    save_time::Float64 = 2.e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "minmod"
    spatial_discretizationScheme_U = "minmod"
    spatial_discretizationScheme_T = "minmod"
    spatial_discretizationScheme_Y = "minmod"
    #spatial_discretizationScheme_Y = "minmod"
    #spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 2440.0)
    push!(γ, 1.648)
    push!(b, 0.0)
    push!(q, 0.0)


    function initial_p(x::Float64,y::Float64)
        if 0.0 < x < 0.05
            return 1.01325e5
        elseif 0.05 <= x < 0.15
            return 1.59060e5
        else
            return 1.01325e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if 0.0 < x < 0.05
            return 0.0
        elseif 0.05 <= x < 0.15
            return 125.65
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64) = 0.0
    function initial_T(x::Float64,y::Float64)
        if 0.0 < x < 0.05
            return 351.82
        elseif 0.05 <= x < 0.15
            return 402.67
        else
            return 351.82
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if x < 0.15
            return 1.0
        else
            return 0.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 125.65
    left_V_BCtype = "fixedValue"
    left_V_BCValue = 0.0
    left_T_BCtype = "fixedValue"
    left_T_BCValue = 402.67
    left_Y_BCtype = "fixedValue"
    left_Y_BCValue = 1.0

    right_P_BCtype = "fixedValue"
    right_P_BCValue = 1.e5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, left_V_BCValue, left_T_BCValue, left_Y_BCValue,
    right_P_BCValue, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end

=#






function oneD_test4_air_water_interface()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-8, 5.e-8, 5.e-8]
    Δt_iters::Vector{UInt32} = [10, 50, 1000000]
    pseudoMaxIter::Vector{UInt32} = [1, 1, 1]
    
    time_end::Float64 = 2.78e-4
    save_time::Float64 = 2.78e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "upwind"
    spatial_discretizationScheme_U = "minmod"
    spatial_discretizationScheme_T = "upwind"
    #spatial_discretizationScheme_Y = "upwind"
    #spatial_discretizationScheme_Y = "minmod"
    spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 621780000.0)
    push!(cᵥ, 3610.0)
    push!(γ, 1.19)
    push!(b, 6.7212e-4)
    push!(q, -1177788.0)


    function initial_p(x::Float64,y::Float64)
        if x < 0.25
            return 1.165e7
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.25
            return 2869.3
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64) = 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.25
            return 6137.3
        else
            return 300.0
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
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 2869.3
    left_V_BCtype = "fixedValue"
    left_V_BCValue = 0.0
    left_T_BCtype = "fixedValue"
    left_T_BCValue = 6137.3
    left_Y_BCtype = "fixedValue"
    left_Y_BCValue = 1.0

    right_P_BCtype = "fixedValue"
    right_P_BCValue = 1.e5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, left_V_BCValue, left_T_BCValue, left_Y_BCValue,
    right_P_BCValue, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end







function oneD_test5_air_water_shock_tube()

    Nx::UInt32 = 800
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 2.0
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.01


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-8, 5.e-8, 1.e-7, 1.e-6]
    Δt_iters::Vector{UInt32} = [10, 20, 100, 1000000]
    pseudoMaxIter::Vector{UInt32} = [3, 3, 3, 3]
    
    time_end::Float64 = 8.e-4
    save_time::Float64 = 8.e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "upwind"
    spatial_discretizationScheme_U = "minmod"
    spatial_discretizationScheme_T = "upwind"
    spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "minmod"
    #spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 621780000.0)
    push!(cᵥ, 3610.0)
    push!(γ, 1.19)
    push!(b, 6.7212e-4)
    push!(q, -1177788.0)


    function initial_p(x::Float64,y::Float64)
        if x < 0.5
            return 1.e9
        else
            return 1.e4
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.5
            return 0.0
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64) = 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.5
            return 300.0
        else
            return 300.0
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if x < 0.5
            return 1.0 - 1.e-5
        else
            return 1.e-5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
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











function oneD_test6_air_helium_shock_bubble()

    Nx::UInt32 = 500
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-5, 1.e-5]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [5, 5]
    
    time_end::Float64 = 6.5e-4
    save_time::Float64 = 6.5e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "minmod"
    spatial_discretizationScheme_U = "minmod"
    spatial_discretizationScheme_T = "minmod"
    spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "minmod"
    #spatial_discretizationScheme_Y = "mstacs"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

    push!(p∞, 0.0)
    push!(cᵥ, 3115.6)
    push!(γ, 1.667)
    push!(b, 0.0)
    push!(q, 0.0)


    function initial_p(x::Float64,y::Float64)
        if x < 0.3
            return 1.245e5
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 0.3
            return 55.33
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64) = 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.3
            return 319.48
        else
            return 300.0
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if 0.5 < x < 0.7
            return 0.0
        else
            return 1.0
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 55.33
    left_V_BCtype = "fixedValue"
    left_V_BCValue = 0.0
    left_T_BCtype = "fixedValue"
    left_T_BCValue = 319.48
    left_Y_BCtype = "fixedValue"
    left_Y_BCValue = 1.0

    right_P_BCtype = "fixedValue"
    right_P_BCValue = 1.e5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, left_V_BCValue, left_T_BCValue, left_Y_BCValue,
    right_P_BCValue, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end







function oneD_test7_air_bubble_in_water()

    Nx::UInt32 = 500
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 2.0
    Ly::Float64 = 0.1
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.1


    Δt::Float64 = 1000.0
    Δt_steps::Vector{Float64} = [1.e-7, 1.e-6]
    Δt_iters::Vector{UInt32} = [10, 1000000]
    pseudoMaxIter::Vector{UInt32} = [1, 1]
    
    time_end::Float64 = 6.5e-4
    save_time::Float64 = 4.0e-4
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    temporal_discretizationScheme = "2nd"
    spatial_discretizationScheme_p = "upwind"
    spatial_discretizationScheme_U = "upwind"
    spatial_discretizationScheme_T = "upwind"
    #spatial_discretizationScheme_Y = "upwind"
    #spatial_discretizationScheme_Y = "minmod"
    spatial_discretizationScheme_Y = "upwind"
    #spatial_discretizationScheme_Y = "superbee"

    gravity = [0.0 0.0 0.0]
    

    push!(p∞, 0.0)
    push!(cᵥ, 720.0)
    push!(γ, 1.4)
    push!(b, 0.0)
    push!(q, 0.0)

#=

    push!(p∞, 621780000.0)
    push!(cᵥ, 3610.0)
    push!(γ, 1.19)
    push!(b, 6.7212e-4)
    push!(q, -1177788.0)
=#

    push!(p∞, 1.804e9)
    push!(cᵥ, 1935.5)
    push!(γ, 4.1)
    push!(b, 0.0)
    push!(q, 0.0)

#=
    push!(p∞, 6.e8)
    push!(cᵥ, 1569.0)
    push!(γ, 4.4)
    push!(b, 0.0)
    push!(q, 0.0)
=#

    function initial_p(x::Float64,y::Float64)
        if x < 1.1
            return 1.487e8
        else
            return 1.e5
        end
    end
    function initial_u(x::Float64,y::Float64)
        if x < 1.1
            return 100.45
        else
            return 0.0
        end
    end
    initial_v(x::Float64,y::Float64) = 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 1.1
            return 302.61
        else
            return 300.0
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if 1.3 < x < 1.5
            return 1.0 - 1.e-5
        else
            return 1.e-5
        end
    end


    left_P_BCtype = "zeroGradient"
    left_U_BCtype = "fixedValue"
    left_U_BCValue = 100.45
    left_V_BCtype = "fixedValue"
    left_V_BCValue = 0.0
    left_T_BCtype = "fixedValue"
    left_T_BCValue = 302.61
    left_Y_BCtype = "fixedValue"
    left_Y_BCValue = 0.0

    right_P_BCtype = "fixedValue"
    right_P_BCValue = 1.e5
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
    temporal_discretizationScheme, 
    spatial_discretizationScheme_p,
    spatial_discretizationScheme_U,
    spatial_discretizationScheme_T,
    spatial_discretizationScheme_Y,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_P_BCtype, left_U_BCtype, left_V_BCtype, left_T_BCtype, left_Y_BCtype,
    right_P_BCtype, right_U_BCtype, right_V_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_P_BCtype, bottom_U_BCtype, bottom_V_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_P_BCtype, top_U_BCtype, top_V_BCtype, top_T_BCtype, top_Y_BCtype,
    nothing, left_U_BCValue, left_V_BCValue, left_T_BCValue, left_Y_BCValue,
    right_P_BCValue, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end


