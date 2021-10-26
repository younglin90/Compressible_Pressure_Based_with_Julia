


function coupled_Ap_boundary!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    bc,
    p_BCtype, p_BCValue, 
    u_BCtype, u_BCValue, 
    v_BCtype, v_BCValue, 
    T_BCtype, T_BCValue, 
    Y_BCtype, Y_BCValue,
    Ap
    )

    for face in bc

        ρₙ = cells[face.owner].var[👉.ρ]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        Tₙ = cells[face.owner].var[👉.T]

        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        uₙ = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        vₙ = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        wₙ = 0.0#cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pₙ = cells[face.owner].var[👉.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pₙ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pₙ = p_BCValue(👉.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uₙ = cells[face.owner].var[👉.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uₙ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uₙ = uₙ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uₙ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uₙ = u_BCValue(👉.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vₙ = cells[face.owner].var[👉.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vₙ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vₙ = vₙ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vₙ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vₙ = v_BCValue(👉.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tₙ = cells[face.owner].var[👉.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tₙ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tₙ = T_BCValue(👉.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Y₁ₙ = cells[face.owner].var[👉.Y₁]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue(👉.time)
        end
        
        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        ρₙ, Hₜₙ, cₙ = faceEOS!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)
        

        flux = ρₙ * Uₙ * ΔS
        Ap[face.owner] += flux / cells[face.owner].Ω


    end
 


end



























function push_A_conv_diff!(
    A_rows::Array{Int64},
    A_cols::Array{Int64},
    A_vals::Array{Float64},
    AiL::Int64, iL::Int64, jL::Int64,
    AiR::Int64, iR::Int64, jR::Int64,
    convfluxₗ::Float64, difffluxₗ::Float64, 
    convfluxᵣ::Float64, difffluxᵣ::Float64
)
    A_vals[AiL] += ( convfluxₗ + difffluxₗ )
    push!(A_rows, iL)
    push!(A_cols, jL)
    push!(A_vals, convfluxᵣ + difffluxᵣ)
    
    A_vals[AiR] -= ( convfluxᵣ + difffluxᵣ )
    push!(A_rows, iR)
    push!(A_cols, jR)
    push!(A_vals, -( convfluxₗ + difffluxₗ ))
end


function coupled!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face}
    )

    
    B_n = 5
    A_n = B_n * B_n

    A_rows = zeros(Int64, length(cells)*A_n)
    A_cols = zeros(Int64, length(cells)*A_n)
    A_vals = zeros(Float64, length(cells)*A_n)
    B = zeros(Float64, length(cells)*B_n)
    
    diagon = 1

    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        Ω = cell.Ω
        Δt = 👉.Δt
        u = cell.var[👉.u]
        v = cell.var[👉.v]
        ρ = cell.var[👉.ρ]
        Hₜ = cell.var[👉.Hₜ]
        p = cell.var[👉.p]
        Y₁ = cell.var[👉.Y₁]
        ∂Hₜ∂p = cell.var[👉.∂Hₜ∂p]
        ∂Hₜ∂T = cell.var[👉.∂Hₜ∂T]
        ∂ρ∂p = cell.var[👉.∂ρ∂p]
        ∂ρ∂T = cell.var[👉.∂ρ∂T]
        ∂ρ∂Y₁ = cell.var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ = cell.var[👉.∂Hₜ∂Y₁]
        ρⁿ = cell.var[👉.ρⁿ]
        uⁿ = cell.var[👉.uⁿ]
        vⁿ = cell.var[👉.vⁿ]
        Hₜⁿ = cell.var[👉.Hₜⁿ]
        pⁿ = cell.var[👉.pⁿ]
        Y₁ⁿ = cell.var[👉.Y₁ⁿ]
        ρⁿ⁻¹ = cell.var[👉.ρⁿ⁻¹]
        uⁿ⁻¹ = cell.var[👉.uⁿ⁻¹]
        vⁿ⁻¹ = cell.var[👉.vⁿ⁻¹]
        Hₜⁿ⁻¹ = cell.var[👉.Hₜⁿ⁻¹]
        pⁿ⁻¹ = cell.var[👉.pⁿ⁻¹]
        Y₁ⁿ⁻¹ = cell.var[👉.Y₁ⁿ⁻¹]

        c_Δt = 1.0
        if 👉.temporal_discretizationScheme == "1st"  
            c_Δt = 1.0
        elseif 👉.temporal_discretizationScheme == "2nd"
            c_Δt = 1.5
        end
        #println(pⁿ,uⁿ,vⁿ,Hₜⁿ)
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = c_Δt * ( ∂ρ∂p*Ω/Δt )
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = c_Δt * ( ∂ρ∂T*Ω/Δt )
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 5
        A_vals[i] = c_Δt * ( ∂ρ∂Y₁*Ω/Δt )

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = c_Δt * ( ∂ρ∂p*Ω/Δt * u ) - ∂ρ∂p*👉.gravity[1]*Ω

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = c_Δt * ( ρ*Ω/Δt )
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = c_Δt * ( ∂ρ∂T*Ω/Δt * u ) - ∂ρ∂T*👉.gravity[1]*Ω
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 5
        A_vals[i] = c_Δt * ( ∂ρ∂Y₁*Ω/Δt * u ) - ∂ρ∂Y₁*👉.gravity[1]*Ω


        # y-momentum
        #g = -9.8
        #g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = c_Δt * ( ∂ρ∂p*Ω/Δt * v ) - ∂ρ∂p*👉.gravity[2]*Ω
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = c_Δt * ( ρ*Ω/Δt )

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = c_Δt * ( ∂ρ∂T*Ω/Δt * v ) - ∂ρ∂T*👉.gravity[2]*Ω
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 5
        A_vals[i] = c_Δt * ( ∂ρ∂Y₁*Ω/Δt * v ) - ∂ρ∂Y₁*👉.gravity[2]*Ω




        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = c_Δt * ( ∂ρ∂p*Ω/Δt * Hₜ + ∂Hₜ∂p*Ω/Δt * ρ - Ω/Δt )

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = c_Δt * ( u*Ω/Δt * ρ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = c_Δt * ( v*Ω/Δt * ρ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = c_Δt * ( ∂ρ∂T*Ω/Δt * Hₜ + ∂Hₜ∂T*Ω/Δt * ρ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 5
        A_vals[i] = c_Δt * ( ∂ρ∂Y₁*Ω/Δt * Hₜ + ∂Hₜ∂Y₁*Ω/Δt * ρ )



        # mass fraction
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 1
        A_vals[i] = c_Δt * ( ∂ρ∂p*Ω/Δt * Y₁ )

        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 4
        A_vals[i] = c_Δt * ( ∂ρ∂T*Ω/Δt * Y₁ )
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 5
        A_vals[i] = c_Δt * ( ∂ρ∂Y₁*Ω/Δt * Y₁ + Ω/Δt * ρ )
        
        

        # B
        if 👉.temporal_discretizationScheme == "1st"
            B[ijStart + 1] = -(ρ - ρⁿ)*Ω/Δt
            B[ijStart + 2] = -(ρ*u - ρⁿ*uⁿ)*Ω/Δt + ρ*👉.gravity[1]*Ω 
            B[ijStart + 3] = -(ρ*v - ρⁿ*vⁿ)*Ω/Δt + ρ*👉.gravity[2]*Ω 
            B[ijStart + 4] = -(ρ*Hₜ - ρⁿ*Hₜⁿ)*Ω/Δt + (p - pⁿ)*Ω/Δt
            B[ijStart + 5] = -(ρ*Y₁ - ρⁿ*Y₁ⁿ)*Ω/Δt
        elseif 👉.temporal_discretizationScheme == "2nd"
            B[ijStart + 1] = -(1.5*ρ - 2.0*ρⁿ + 0.5*ρⁿ⁻¹)*Ω/Δt
            B[ijStart + 2] = -(1.5*ρ*u - 2.0*ρⁿ*uⁿ + 0.5*ρⁿ⁻¹*uⁿ⁻¹)*Ω/Δt + ρ*👉.gravity[1]*Ω 
            B[ijStart + 3] = -(1.5*ρ*v - 2.0*ρⁿ*vⁿ + 0.5*ρⁿ⁻¹*vⁿ⁻¹)*Ω/Δt + ρ*👉.gravity[2]*Ω 
            B[ijStart + 4] = -(1.5*ρ*Hₜ - 2.0*ρⁿ*Hₜⁿ + 0.5*ρⁿ⁻¹*Hₜⁿ⁻¹)*Ω/Δt + (1.5*p - 2.0*pⁿ + 0.5*pⁿ⁻¹)*Ω/Δt
            B[ijStart + 5] = -(1.5*ρ*Y₁ - 2.0*ρⁿ*Y₁ⁿ + 0.5*ρⁿ⁻¹*Y₁ⁿ⁻¹)*Ω/Δt
        end



        diagon += 1


    end

    
    ∂Δp∂x0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pₙ = 0.5 * (cells[face.owner].var[👉.p] + cells[face.neighbour].var[👉.p])
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.neighbour, 1] -= pₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 2] -= pₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 3] -= pₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end

    for face in faces_boundary
        pₙ = cells[face.owner].var[👉.p]
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end






    
    Ap = zeros(Float64, length(cells))
    for face in faces_internal


        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        #uₗ = cells[face.owner].var[👉.u]
        #uᵣ = cells[face.neighbour].var[👉.u]
        #vₗ = cells[face.owner].var[👉.v]
        #vᵣ = cells[face.neighbour].var[👉.v]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d̂ = 👉.Δt / ρˢ
        
        Wₗ = 0.0
        Wᵣ = 0.0
        #if 👉.spatial_discretizationScheme == "upwind"
            Wₗ = 0.5 * (1.0 + sign(Uₙ))
            Wᵣ = 1.0 - Wₗ
        #elseif 👉.spatial_discretizationScheme == "central"
        #    Wₗ = 0.5
        #    Wᵣ = 1.0 - Wₗ
        #end

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        # Rhie-Chow
        #=
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR
        =#

        flux = ρₙ * Uₙ * face.ΔS
        Ap[face.owner] += Wₗ * flux / cells[face.owner].Ω
        Ap[face.neighbour] -= Wᵣ * flux / cells[face.neighbour].Ω
    end

    coupled_Ap_boundary!(
    👉,cells,faces,
    faces_boundary_top, 
    👉.top_p_BCtype, 👉.top_p_BCValue, 
    👉.top_u_BCtype, 👉.top_u_BCValue, 
    👉.top_v_BCtype, 👉.top_v_BCValue, 
    👉.top_T_BCtype, 👉.top_T_BCValue, 
    👉.top_Y_BCtype, 👉.top_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(👉,cells,faces,
    faces_boundary_bottom, 
    👉.bottom_p_BCtype, 👉.bottom_p_BCValue, 
    👉.bottom_u_BCtype, 👉.bottom_u_BCValue, 
    👉.bottom_v_BCtype, 👉.bottom_v_BCValue, 
    👉.bottom_T_BCtype, 👉.bottom_T_BCValue, 
    👉.bottom_Y_BCtype, 👉.bottom_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(👉,cells,faces,
    faces_boundary_left, 
    👉.left_p_BCtype, 👉.left_p_BCValue, 
    👉.left_u_BCtype, 👉.left_u_BCValue, 
    👉.left_v_BCtype, 👉.left_v_BCValue, 
    👉.left_T_BCtype, 👉.left_T_BCValue, 
    👉.left_Y_BCtype, 👉.left_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(👉,cells,faces,
    faces_boundary_right, 
    👉.right_p_BCtype, 👉.right_p_BCValue, 
    👉.right_u_BCtype, 👉.right_u_BCValue, 
    👉.right_v_BCtype, 👉.right_v_BCValue, 
    👉.right_T_BCtype, 👉.right_T_BCValue, 
    👉.right_Y_BCtype, 👉.right_Y_BCValue,
    Ap)
    






    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

       # ρₗ = cells[face.owner].var[👉.ρ]
       # ρᵣ = cells[face.neighbour].var[👉.ρ]

        ρₗ = face.varₗ[👉.ρ]
        ρᵣ = face.varᵣ[👉.ρ]

        pO = cells[face.owner].var[👉.p]
        pN = cells[face.neighbour].var[👉.p]
        pₗ = face.varₗ[👉.p]
        pᵣ = face.varᵣ[👉.p]

       # uₗ = cells[face.owner].var[👉.u]
       # uᵣ = cells[face.neighbour].var[👉.u]
       # vₗ = cells[face.owner].var[👉.v]
       # vᵣ = cells[face.neighbour].var[👉.v]

        uₗ = face.varₗ[👉.u]
        uᵣ = face.varᵣ[👉.u]
        vₗ = face.varₗ[👉.v]
        vᵣ = face.varᵣ[👉.v]

        wₗ = 0.0#cells[face.owner].var[👉.w]
        wᵣ = 0.0#cells[face.neighbour].var[👉.w]
       #= Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
        #μₗ = cells[face.owner].var[👉.μ]
        #μᵣ = cells[face.neighbour].var[👉.μ]
        ∂ρ∂pₗ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂pᵣ = cells[face.neighbour].var[👉.∂ρ∂p]
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = cells[face.neighbour].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]
        Y₁ₗ = cells[face.owner].var[👉.Y₁]
        Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]
        ∂ρ∂Y₁ₗ = cells[face.owner].var[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ₗ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂Hₜ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂Hₜ∂Y₁]
        =#
        
        Hₜₗ = face.varₗ[👉.Hₜ]
        Hₜᵣ = face.varᵣ[👉.Hₜ]
        ∂ρ∂pₗ = face.varₗ[👉.∂ρ∂p]
        ∂ρ∂pᵣ = face.varᵣ[👉.∂ρ∂p]
        ∂ρ∂Tₗ = face.varₗ[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = face.varᵣ[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = face.varₗ[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = face.varᵣ[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = face.varₗ[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = face.varᵣ[👉.∂Hₜ∂T]
        Y₁ₗ = face.varₗ[👉.Y₁]
        Y₁ᵣ = face.varᵣ[👉.Y₁]
        ∂ρ∂Y₁ₗ = face.varₗ[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = face.varᵣ[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ₗ = face.varₗ[👉.∂Hₜ∂Y₁]
        ∂Hₜ∂Y₁ᵣ = face.varᵣ[👉.∂Hₜ∂Y₁]
        cₗ = face.varₗ[👉.c]
        cᵣ = face.varᵣ[👉.c]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        #d = 0.5 * (1.0 / (Ap[face.owner]) + 1.0 / (Ap[face.neighbour]) )
        d̂ = 👉.Δt / ρˢ
        #d̂ = d / (2.0 + ρˢ / 👉.Δt * d)
        #if d>1.e9
        #    d̂ = 👉.Δt / ρˢ
        #end
        
        # Rhie-Chow
        Uₙ_RC = 0.0
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ_RC += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ_RC -= d̂ * (pN-pO) / ΔLR

        RCdiffΔp = d̂ / ΔLR

        # before step
        ρₗⁿ = cells[face.owner].var[👉.ρⁿ]
        ρᵣⁿ = cells[face.neighbour].var[👉.ρⁿ]
        ρˢⁿ = 1.0 / (0.5/ρₗⁿ + 0.5/ρᵣⁿ)
        Uₙₗⁿ = cells[face.owner].var[👉.uⁿ] * face.n̂[1] + cells[face.owner].var[👉.vⁿ] * face.n̂[2]
        Uₙᵣⁿ = cells[face.neighbour].var[👉.uⁿ] * face.n̂[1] + cells[face.neighbour].var[👉.vⁿ] * face.n̂[2]
        #Uₙ_RC += d̂ * ρˢⁿ / 👉.Δt * ( face.Uₙⁿ - 0.5 * (Uₙₗⁿ + Uₙᵣⁿ) )
        #Uₙ += ( face.Uₙⁿ - 0.5 * (Uₙₗⁿ + Uₙᵣⁿ) )

        # YYL riemann
        c̄ = 0.5*(cₗ + cᵣ)
        Mₗ = Uₙₗ/c̄
        Mᵣ = Uₙᵣ/c̄
        # calculate M+ and P+ for left state
        Mₗ⁺ = M_func(Mₗ,1.0,0.125)
        p⁺ = pre_func(Mₗ,1.0,0.1875)
        # calculate M- and P- for left state
        Mᵣ⁻ = M_func(Mᵣ,-1.0,0.125)
        p⁻ = pre_func(Mᵣ,-1.0,0.1875)
        KLR = sqrt(0.5*(uₗ^2+vₗ^2+wₗ^2+uᵣ^2+vᵣ^2+wᵣ^2))
        Mdash = min(1.0,KLR/c̄)

        CdiffΔp = 0.5*(1.0-Mdash)^2/c̄ /(0.5*(ρₗ+ρᵣ))
        WUₗ = Mₗ⁺/Mₗ
        if abs(Mₗ) <= 0.00001
            WUₗ = 0.5
        end
        WUᵣ = Mᵣ⁻/Mᵣ
        if abs(Mᵣ) <= 0.00001
            WUᵣ = 0.5
        end
        WUₗ = 0.5
        WUᵣ = 0.5

	    Uₙ = WUₗ*Uₙₗ + WUᵣ*Uₙᵣ + Uₙ_RC #- CdiffΔp * (pᵣ-pₗ)
	    #Uₙ = WUₗ*Uₙₗ + WUᵣ*Uₙᵣ - CdiffΔp * (pᵣ-pₗ)

        Wpₗ = p⁺
        Wpᵣ = p⁻

        diffΔp = RCdiffΔp #+ CdiffΔp
        #diffΔp = CdiffΔp

        #--------------------
        # SAVE
        face.Uₙ = Uₙ
        #--------------------




        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Hₜₙ = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ
        Y₁ₙ = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ

        pₙ = Wpₗ*pₗ + Wpᵣ*pᵣ

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)


        #------------------------
        # continuity
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * diffΔp * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * diffΔp * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * diffΔp * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * diffΔp * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        
        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[1] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[2] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS ))
        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Uₙ * ΔS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + Wpₗ * face.n̂[1] * ΔS + ρₙ * diffΔp * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + Wpᵣ * face.n̂[1] * ΔS - ρₙ * diffΔp * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + Wpᵣ * face.n̂[1] * ΔS - ρₙ * diffΔp * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + Wpₗ * face.n̂[1] * ΔS + ρₙ * diffΔp * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[2] * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS ))

        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * uₙ * Uₙ * ΔS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] +=  ( Wₗ * ∂ρ∂pₗ* vₙ * Uₙ * ΔS + Wpₗ * face.n̂[2] * ΔS + ρₙ * diffΔp * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + Wpᵣ * face.n̂[2] * ΔS - ρₙ * diffΔp * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ* vₙ * Uₙ * ΔS + Wpᵣ * face.n̂[2] * ΔS - ρₙ * diffΔp * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + Wpₗ * face.n̂[2] * ΔS + ρₙ * diffΔp * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[1] * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS ))

        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -(  Wₗ * ∂ρ∂Y₁ₗ * vₙ * Uₙ *ΔS ))


        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * diffΔp * Hₜₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * diffΔp * Hₜₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * diffΔp * Hₜₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * diffΔp * Hₜₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS ))

        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Y₁ᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Y₁ₗ * Uₙ * ΔS ))
        

        

        #------------------------
        # mass fraction
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * diffΔp * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Y₁ₙ * Uₙ * ΔS - ρₙ * diffΔp * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Y₁ₙ * Uₙ * ΔS - ρₙ * diffΔp * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * diffΔp * Y₁ₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[1] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[1] * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[1] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[1] * Y₁ₙ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * WUₗ * face.n̂[2] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * WUᵣ * face.n̂[2] * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * WUᵣ * face.n̂[2] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * WUₗ * face.n̂[2] * Y₁ₙ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Y₁ₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Y₁ₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Y₁ₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Y₁ₙ * Uₙ * ΔS ))

        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wₗ * Uₙ * ΔS ))
        

        # ----------------------------

        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartᵣ + 1] += ( ρₙ * Uₙ * ΔS )

        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )

        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )

        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartᵣ + 4] += ( ρₙ * Hₜₙ * Uₙ * ΔS )

        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        B[ijStartᵣ + 5] += ( ρₙ * Y₁ₙ * Uₙ * ΔS )




    end


    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    #bc_wall = []
    #append!( bc_wall, faces_boundary_top )
    #append!( bc_wall, faces_boundary_bottom )
    #append!( bc_wall, faces_boundary_left )
    #append!( bc_wall, faces_boundary_right )


    coupled_boundary!(
        👉,cells,faces,
        faces_boundary_top, 
        👉.top_p_BCtype, 👉.top_p_BCValue, 
        👉.top_u_BCtype, 👉.top_u_BCValue, 
        👉.top_v_BCtype, 👉.top_v_BCValue, 
        👉.top_T_BCtype, 👉.top_T_BCValue, 
        👉.top_Y_BCtype, 👉.top_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(👉,cells,faces,
        faces_boundary_bottom, 
        👉.bottom_p_BCtype, 👉.bottom_p_BCValue, 
        👉.bottom_u_BCtype, 👉.bottom_u_BCValue, 
        👉.bottom_v_BCtype, 👉.bottom_v_BCValue, 
        👉.bottom_T_BCtype, 👉.bottom_T_BCValue, 
        👉.bottom_Y_BCtype, 👉.bottom_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(👉,cells,faces,
        faces_boundary_left, 
        👉.left_p_BCtype, 👉.left_p_BCValue, 
        👉.left_u_BCtype, 👉.left_u_BCValue, 
        👉.left_v_BCtype, 👉.left_v_BCValue, 
        👉.left_T_BCtype, 👉.left_T_BCValue, 
        👉.left_Y_BCtype, 👉.left_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(👉,cells,faces,
        faces_boundary_right, 
        👉.right_p_BCtype, 👉.right_p_BCValue, 
        👉.right_u_BCtype, 👉.right_u_BCValue, 
        👉.right_v_BCtype, 👉.right_v_BCValue, 
        👉.right_T_BCtype, 👉.right_T_BCValue, 
        👉.right_Y_BCtype, 👉.right_Y_BCValue,
        B_n, A_n, A_vals, B)


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ΔQ = solve(ps, A, B)

    #ΔQ = A\B

    #ml = ruge_stuben(A)
    #ΔQ = solve(ml, A)
    #P = aspreconditioner(ml)
    #ΔQ = bicgstabl(A, B, Pl = P)
    #ΔQ = gmres(A, B)

    #ΔQ = A\B

   # error()
   # exit()

   relax_p = 1.0
   relax_U = 1.0
   relax_T = 1.0
   relax_Y = 1.0


   diagon = 1
   maximum_p = -1.e12
   maximum_U = -1.e12
   maximum_T = -1.e12
   maximum_Y = -1.e12
   norm_p = 0.0
   norm_U = 0.0
   norm_T = 0.0
   norm_total = 0.0
   for cell in cells

       ijStart = B_n*(diagon-1)
       Astart = A_n*(diagon-1)
       i = Astart

       cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
       cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
       cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
       cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]
       cell.var[👉.Y₁] += relax_Y * ΔQ[ijStart + 5]

       cell.var[👉.p] = max(cell.var[👉.p],1.e-200)
       cell.var[👉.T] = max(cell.var[👉.T],1.e-200)
       
       norm_total += cell.var[👉.p]^2
       norm_total += cell.var[👉.u]^2
       norm_total += cell.var[👉.v]^2
       norm_total += cell.var[👉.T]^2
       norm_total += cell.var[👉.Y₁]^2

       diagon += 1
   end

   norm_total = sqrt(norm_total)

   #sleep(1000.0)


   return norm(ΔQ)/length(cells)/5
  


end



function coupled_boundary!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    bc,
    p_BCtype, p_BCValue, 
    u_BCtype, u_BCValue, 
    v_BCtype, v_BCValue, 
    T_BCtype, T_BCValue, 
    Y_BCtype, Y_BCValue,
    B_n, A_n, A_vals, B
    )

    for face in bc
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        Tₙ = cells[face.owner].var[👉.T]

        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        uₙ = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        vₙ = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        wₙ = 0.0#cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        id = []
        push!(id,i)
        push!(id,i+5)
        push!(id,i+10)
        push!(id,i+15)
        push!(id,i+20)

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pₙ = cells[face.owner].var[👉.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pₙ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pₙ = p_BCValue(👉.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uₙ = cells[face.owner].var[👉.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uₙ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uₙ = uₙ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uₙ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uₙ = u_BCValue(👉.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vₙ = cells[face.owner].var[👉.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vₙ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vₙ = vₙ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vₙ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vₙ = v_BCValue(👉.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tₙ = cells[face.owner].var[👉.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tₙ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tₙ = T_BCValue(👉.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Y₁ₙ = cells[face.owner].var[👉.Y₁]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue(👉.time)
        end
        
        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        ρₙ, Hₜₙ, cₙ = faceEOS!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)
        
        # continuity
        i += 1
        A_vals[i] += coeff_p * (∂ρ∂pₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_u * (ρₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += coeff_v * (ρₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += coeff_T * (∂ρ∂Tₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_Y * (∂ρ∂Y₁ₙ * Uₙ * ΔS)

        
        # x-momentum
        i += 1
        A_vals[i] += coeff_p * (∂ρ∂pₙ * uₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_u * (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += coeff_v * (ρₙ * uₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += coeff_T * (∂ρ∂Tₙ * uₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_Y * (∂ρ∂Y₁ₙ * uₙ * Uₙ * ΔS)

        
        # y-momentum
        i += 1
        A_vals[i] += coeff_p * (∂ρ∂pₙ * vₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_u * (ρₙ * vₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += coeff_v * (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += coeff_T * (∂ρ∂Tₙ * vₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_Y * (∂ρ∂Y₁ₙ * vₙ * Uₙ * ΔS)


        # energy
        i += 1
        A_vals[i] += coeff_p * (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS)
        i += 1
        A_vals[i] += coeff_u * (ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS)
        i += 1
        A_vals[i] += coeff_v * (ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS)
        i += 1
        A_vals[i] += coeff_T * (∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS)
        i += 1
        A_vals[i] += coeff_Y * (∂ρ∂Y₁ₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Y₁ₙ * ΔS)


        # massfraction
        i += 1
        A_vals[i] += coeff_p * (∂ρ∂pₙ * Uₙ * Y₁ₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += coeff_u * (ρₙ * face.n̂[1] * Y₁ₙ * ΔS)
        i += 1
        A_vals[i] += coeff_v * (ρₙ * face.n̂[2] * Y₁ₙ * ΔS)
        i += 1
        A_vals[i] += coeff_T * (∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS)
        i += 1
        A_vals[i] += coeff_Y * (∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS)


        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        

    end
 


end












function coupled_Ap_boundary!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    bc,
    p_BCtype, p_BCValue, 
    u_BCtype, u_BCValue, 
    v_BCtype, v_BCValue, 
    T_BCtype, T_BCValue, 
    Y_BCtype, Y_BCValue,
    Ap
    )

    for face in bc

        ρₙ = cells[face.owner].var[👉.ρ]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        Tₙ = cells[face.owner].var[👉.T]

        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        uₙ = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        vₙ = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        wₙ = 0.0#cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pₙ = cells[face.owner].var[👉.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pₙ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pₙ = p_BCValue(👉.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uₙ = cells[face.owner].var[👉.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uₙ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uₙ = uₙ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uₙ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uₙ = u_BCValue(👉.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vₙ = cells[face.owner].var[👉.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vₙ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vₙ = vₙ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vₙ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vₙ = v_BCValue(👉.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tₙ = cells[face.owner].var[👉.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tₙ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tₙ = T_BCValue(👉.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Y₁ₙ = cells[face.owner].var[👉.Y₁]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Y₁ₙ = Y_BCValue(👉.time)
        end
        
        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        ρₙ, Hₜₙ, cₙ = faceEOS!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)
        

        flux = ρₙ * Uₙ * ΔS
        Ap[face.owner] += flux / cells[face.owner].Ω


    end
 


end


function M_func(M::Float64, op::Float64, α::Float64)
    mu = 0.0
	if abs(M) > 1.0 
		mu = 0.5*(M + op*abs(M))
	else
		mu = op*0.25*(M + op)^2.0 + op*α*(M*M-1.0)^2.0
    end
	
	return mu
end

function pre_func(M::Float64, op::Float64, α::Float64)
    mu = 0.0
	if abs(M) > 1.0
		mu = 0.5*(1.0 + op*sign(M) )
	else
		mu = 0.25*(M + op)^2.0*(2.0-op*M) + op*α*M*(M*M-1.0)^2.0
    end
	
	return mu;
end


