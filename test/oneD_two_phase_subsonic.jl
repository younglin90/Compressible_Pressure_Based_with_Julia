

function oneD_interface_advection_with_constant_velocity()

    Nx::UInt32 = 100
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    Δt::Float64 = 1.e-4
    pseudoMaxIter::UInt32 = 3

    initial_p(x::Float64,y::Float64) = 1.e5
    initial_u(x::Float64,y::Float64) = 1.0
    initial_v(x::Float64,y::Float64) = 0.0
    initial_T(x::Float64,y::Float64) = 300.0
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
    
    return Nx, Ny, Nz, Lx, Ly, Lz, Δt, pseudoMaxIter,
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




function oneD_helium_bubble_in_air()

    Nx::UInt32 = 100
    Ny::UInt32 = 1
    Nz::UInt32 = 1
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0
    Lz::Float64 = 0.1
    Δt::Float64 = 1.e-4
    pseudoMaxIter::UInt32 = 3

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
    initial_v(x::Float64,y::Float64)= 0.0
    function initial_T(x::Float64,y::Float64)
        if x < 0.3
            return 319.48
        else
            return 300.0
        end
    end
    function initial_Y(x::Float64,y::Float64)
        if 0.5 < x < 0.7
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
    
    return Nx, Ny, Nz, Lx, Ly, Lz, Δt, pseudoMaxIter,
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
