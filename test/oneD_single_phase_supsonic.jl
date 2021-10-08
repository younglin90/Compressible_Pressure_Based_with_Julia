
function sod_shock_test()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    Δt::Float64 = 3.e-4
    pseudoMaxIter::UInt32 = 5

    time_end::Float64 = 0.2
    save_time::Float64 = 0.2
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    discretizationScheme = "upwind"

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
            return 1.0
        else
            return 0.1
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
            return 0.003484
        else
            return 0.002787
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
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing

end





function strong_shock_tube()

    Nx::UInt32 = 1000
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    corantNumber::Float64 = 0.5


    Δt::Float64 = 1.e-6
    Δt_steps::Vector{Float64} = [1.e-6, 5.e-6, 8.e-6]
    Δt_iters::Vector{UInt32} = [50, 100, 1000000]
    pseudoMaxIter::Vector{UInt32} = [7, 8, 9]
    
    time_end::Float64 = 0.012
    save_time::Float64 = 0.012
    save_iteration::UInt32 = 10000

    p∞::Array{Float64} = []
	cᵥ::Array{Float64} = []
    γ::Array{Float64} = []
    b::Array{Float64} = []
	q::Array{Float64} = []

    discretizationScheme = "upwind"

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
            return 1000.0
        else
            return 0.01
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
            return 3.4722222222222222222222222222222
        else
            return 3.4722222222222222222222222222222e-5
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
    discretizationScheme,
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
