


function coupled_Ap_boundary!(
    ð::controls,
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

        Ïâ = cells[face.owner].var[ð.Ï]
        pâ = cells[face.owner].var[ð.p]
        Hââ = cells[face.owner].var[ð.Hâ]
        Yââ = cells[face.owner].var[ð.Yâ]
        Tâ = cells[face.owner].var[ð.T]

        ÎS = face.ÎS

        Uâ = 0.0
        Uâ += cells[face.owner].var[ð.u]*face.nÌ[1]
        Uâ += cells[face.owner].var[ð.v]*face.nÌ[2]
        Uâ += cells[face.owner].var[ð.w]*face.nÌ[3]

        uâ = cells[face.owner].var[ð.u] - Uâ * face.nÌ[1]
        vâ = cells[face.owner].var[ð.v] - Uâ * face.nÌ[2]
        wâ = 0.0#cells[face.owner].var[ð.w] - Uâ * face.nÌ[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pâ = cells[face.owner].var[ð.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pâ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pâ = p_BCValue(ð.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uâ = cells[face.owner].var[ð.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uâ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uâ = uâ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uâ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uâ = u_BCValue(ð.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vâ = cells[face.owner].var[ð.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vâ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vâ = vâ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vâ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vâ = v_BCValue(ð.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tâ = cells[face.owner].var[ð.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tâ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tâ = T_BCValue(ð.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yââ = cells[face.owner].var[ð.Yâ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yââ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yââ = Y_BCValue(ð.time)
        end
        
        Uâ = uâ * face.nÌ[1] + vâ * face.nÌ[2] + wâ * face.nÌ[3]

        Ïâ, Hââ, câ = faceEOS!(ð,pâ,uâ,vâ,wâ,Tâ,Yââ)
        

        flux = Ïâ * Uâ * ÎS
        Ap[face.owner] += flux / cells[face.owner].Î©


    end
 


end



























function push_A_conv_diff!(
    A_rows::Array{Int64},
    A_cols::Array{Int64},
    A_vals::Array{Float64},
    AiL::Int64, iL::Int64, jL::Int64,
    AiR::Int64, iR::Int64, jR::Int64,
    convfluxâ::Float64, difffluxâ::Float64, 
    convfluxáµ£::Float64, difffluxáµ£::Float64
)
    A_vals[AiL] += ( convfluxâ + difffluxâ )
    push!(A_rows, iL)
    push!(A_cols, jL)
    push!(A_vals, convfluxáµ£ + difffluxáµ£)
    
    A_vals[AiR] -= ( convfluxáµ£ + difffluxáµ£ )
    push!(A_rows, iR)
    push!(A_cols, jR)
    push!(A_vals, -( convfluxâ + difffluxâ ))
end


function coupled!(
    ð::controls,
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

        Î© = cell.Î©
        Ît = ð.Ît
        u = cell.var[ð.u]
        v = cell.var[ð.v]
        Ï = cell.var[ð.Ï]
        Hâ = cell.var[ð.Hâ]
        p = cell.var[ð.p]
        Yâ = cell.var[ð.Yâ]
        âHââp = cell.var[ð.âHââp]
        âHââT = cell.var[ð.âHââT]
        âÏâp = cell.var[ð.âÏâp]
        âÏâT = cell.var[ð.âÏâT]
        âÏâYâ = cell.var[ð.âÏâYâ]
        âHââYâ = cell.var[ð.âHââYâ]
        Ïâ¿ = cell.var[ð.Ïâ¿]
        uâ¿ = cell.var[ð.uâ¿]
        vâ¿ = cell.var[ð.vâ¿]
        Hââ¿ = cell.var[ð.Hââ¿]
        pâ¿ = cell.var[ð.pâ¿]
        Yââ¿ = cell.var[ð.Yââ¿]
        Ïâ¿â»Â¹ = cell.var[ð.Ïâ¿â»Â¹]
        uâ¿â»Â¹ = cell.var[ð.uâ¿â»Â¹]
        vâ¿â»Â¹ = cell.var[ð.vâ¿â»Â¹]
        Hââ¿â»Â¹ = cell.var[ð.Hââ¿â»Â¹]
        pâ¿â»Â¹ = cell.var[ð.pâ¿â»Â¹]
        Yââ¿â»Â¹ = cell.var[ð.Yââ¿â»Â¹]

        c_Ît = 1.0
        if ð.temporal_discretizationScheme == "1st"  
            c_Ît = 1.0
        elseif ð.temporal_discretizationScheme == "2nd"
            c_Ît = 1.5
        end
        #println(pâ¿,uâ¿,vâ¿,Hââ¿)
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = c_Ît * ( âÏâp*Î©/Ît )
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = c_Ît * ( âÏâT*Î©/Ît )
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 5
        A_vals[i] = c_Ît * ( âÏâYâ*Î©/Ît )

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = c_Ît * ( âÏâp*Î©/Ît * u ) - âÏâp*ð.gravity[1]*Î©

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = c_Ît * ( Ï*Î©/Ît )
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = c_Ît * ( âÏâT*Î©/Ît * u ) - âÏâT*ð.gravity[1]*Î©
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 5
        A_vals[i] = c_Ît * ( âÏâYâ*Î©/Ît * u ) - âÏâYâ*ð.gravity[1]*Î©


        # y-momentum
        #g = -9.8
        #g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = c_Ît * ( âÏâp*Î©/Ît * v ) - âÏâp*ð.gravity[2]*Î©
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = c_Ît * ( Ï*Î©/Ît )

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = c_Ît * ( âÏâT*Î©/Ît * v ) - âÏâT*ð.gravity[2]*Î©
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 5
        A_vals[i] = c_Ît * ( âÏâYâ*Î©/Ît * v ) - âÏâYâ*ð.gravity[2]*Î©




        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = c_Ît * ( âÏâp*Î©/Ît * Hâ + âHââp*Î©/Ît * Ï - Î©/Ît )

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = c_Ît * ( u*Î©/Ît * Ï )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = c_Ît * ( v*Î©/Ît * Ï )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = c_Ît * ( âÏâT*Î©/Ît * Hâ + âHââT*Î©/Ît * Ï )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 5
        A_vals[i] = c_Ît * ( âÏâYâ*Î©/Ît * Hâ + âHââYâ*Î©/Ît * Ï )



        # mass fraction
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 1
        A_vals[i] = c_Ît * ( âÏâp*Î©/Ît * Yâ )

        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 4
        A_vals[i] = c_Ît * ( âÏâT*Î©/Ît * Yâ )
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 5
        A_vals[i] = c_Ît * ( âÏâYâ*Î©/Ît * Yâ + Î©/Ît * Ï )
        
        

        # B
        if ð.temporal_discretizationScheme == "1st"
            B[ijStart + 1] = -(Ï - Ïâ¿)*Î©/Ît
            B[ijStart + 2] = -(Ï*u - Ïâ¿*uâ¿)*Î©/Ît + Ï*ð.gravity[1]*Î© 
            B[ijStart + 3] = -(Ï*v - Ïâ¿*vâ¿)*Î©/Ît + Ï*ð.gravity[2]*Î© 
            B[ijStart + 4] = -(Ï*Hâ - Ïâ¿*Hââ¿)*Î©/Ît + (p - pâ¿)*Î©/Ît
            B[ijStart + 5] = -(Ï*Yâ - Ïâ¿*Yââ¿)*Î©/Ît
        elseif ð.temporal_discretizationScheme == "2nd"
            B[ijStart + 1] = -(1.5*Ï - 2.0*Ïâ¿ + 0.5*Ïâ¿â»Â¹)*Î©/Ît
            B[ijStart + 2] = -(1.5*Ï*u - 2.0*Ïâ¿*uâ¿ + 0.5*Ïâ¿â»Â¹*uâ¿â»Â¹)*Î©/Ît + Ï*ð.gravity[1]*Î© 
            B[ijStart + 3] = -(1.5*Ï*v - 2.0*Ïâ¿*vâ¿ + 0.5*Ïâ¿â»Â¹*vâ¿â»Â¹)*Î©/Ît + Ï*ð.gravity[2]*Î© 
            B[ijStart + 4] = -(1.5*Ï*Hâ - 2.0*Ïâ¿*Hââ¿ + 0.5*Ïâ¿â»Â¹*Hââ¿â»Â¹)*Î©/Ît + (1.5*p - 2.0*pâ¿ + 0.5*pâ¿â»Â¹)*Î©/Ît
            B[ijStart + 5] = -(1.5*Ï*Yâ - 2.0*Ïâ¿*Yââ¿ + 0.5*Ïâ¿â»Â¹*Yââ¿â»Â¹)*Î©/Ît
        end



        diagon += 1


    end

    
    âÎpâx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pâ = 0.5 * (cells[face.owner].var[ð.p] + cells[face.neighbour].var[ð.p])
        âÎpâx0[face.owner, 1] += pâ * face.nÌ[1] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 2] += pâ * face.nÌ[2] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 3] += pâ * face.nÌ[3] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.neighbour, 1] -= pâ * face.nÌ[1] * face.ÎS / cells[face.neighbour].Î©
        âÎpâx0[face.neighbour, 2] -= pâ * face.nÌ[2] * face.ÎS / cells[face.neighbour].Î©
        âÎpâx0[face.neighbour, 3] -= pâ * face.nÌ[3] * face.ÎS / cells[face.neighbour].Î©
    end

    for face in faces_boundary
        pâ = cells[face.owner].var[ð.p]
        âÎpâx0[face.owner, 1] += pâ * face.nÌ[1] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 2] += pâ * face.nÌ[2] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 3] += pâ * face.nÌ[3] * face.ÎS / cells[face.owner].Î©
    end






    
    Ap = zeros(Float64, length(cells))
    for face in faces_internal


        Ïâ = cells[face.owner].var[ð.Ï]
        Ïáµ£ = cells[face.neighbour].var[ð.Ï]
        pâ = cells[face.owner].var[ð.p]
        páµ£ = cells[face.neighbour].var[ð.p]
        #uâ = cells[face.owner].var[ð.u]
        #uáµ£ = cells[face.neighbour].var[ð.u]
        #vâ = cells[face.owner].var[ð.v]
        #váµ£ = cells[face.neighbour].var[ð.v]
        uâ = cells[face.owner].var[ð.u]
        uáµ£ = cells[face.neighbour].var[ð.u]
        vâ = cells[face.owner].var[ð.v]
        váµ£ = cells[face.neighbour].var[ð.v]

        Uââ = uâ * face.nÌ[1] + vâ * face.nÌ[2]
        Uâáµ£ = uáµ£ * face.nÌ[1] + váµ£ * face.nÌ[2]
        Uâ = 0.5 * (Uââ + Uâáµ£)

        centerâ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centeráµ£ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ÎLR = norm(centeráµ£ - centerâ)

        ÏË¢ = 1.0 / (0.5/Ïâ + 0.5/Ïáµ£)
        dÌ = ð.Ît / ÏË¢
        
        Wâ = 0.0
        Wáµ£ = 0.0
        #if ð.spatial_discretizationScheme == "upwind"
            Wâ = 0.5 * (1.0 + sign(Uâ))
            Wáµ£ = 1.0 - Wâ
        #elseif ð.spatial_discretizationScheme == "central"
        #    Wâ = 0.5
        #    Wáµ£ = 1.0 - Wâ
        #end

        Ïâ = Wâ * Ïâ + Wáµ£ * Ïáµ£
        
        # Rhie-Chow
        #=
        Uâ += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 1] * face.nÌ[1]
        Uâ += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 2] * face.nÌ[2]
        Uâ += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 3] * face.nÌ[3]
        Uâ += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 1] * face.nÌ[1]
        Uâ += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 2] * face.nÌ[2]
        Uâ += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 3] * face.nÌ[3]
        Uâ -= dÌ * (páµ£-pâ) / ÎLR
        =#

        flux = Ïâ * Uâ * face.ÎS
        Ap[face.owner] += Wâ * flux / cells[face.owner].Î©
        Ap[face.neighbour] -= Wáµ£ * flux / cells[face.neighbour].Î©
    end

    coupled_Ap_boundary!(
    ð,cells,faces,
    faces_boundary_top, 
    ð.top_p_BCtype, ð.top_p_BCValue, 
    ð.top_u_BCtype, ð.top_u_BCValue, 
    ð.top_v_BCtype, ð.top_v_BCValue, 
    ð.top_T_BCtype, ð.top_T_BCValue, 
    ð.top_Y_BCtype, ð.top_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(ð,cells,faces,
    faces_boundary_bottom, 
    ð.bottom_p_BCtype, ð.bottom_p_BCValue, 
    ð.bottom_u_BCtype, ð.bottom_u_BCValue, 
    ð.bottom_v_BCtype, ð.bottom_v_BCValue, 
    ð.bottom_T_BCtype, ð.bottom_T_BCValue, 
    ð.bottom_Y_BCtype, ð.bottom_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(ð,cells,faces,
    faces_boundary_left, 
    ð.left_p_BCtype, ð.left_p_BCValue, 
    ð.left_u_BCtype, ð.left_u_BCValue, 
    ð.left_v_BCtype, ð.left_v_BCValue, 
    ð.left_T_BCtype, ð.left_T_BCValue, 
    ð.left_Y_BCtype, ð.left_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(ð,cells,faces,
    faces_boundary_right, 
    ð.right_p_BCtype, ð.right_p_BCValue, 
    ð.right_u_BCtype, ð.right_u_BCValue, 
    ð.right_v_BCtype, ð.right_v_BCValue, 
    ð.right_T_BCtype, ð.right_T_BCValue, 
    ð.right_Y_BCtype, ð.right_Y_BCValue,
    Ap)
    






    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartâ = B_n*(face.owner-1)
        ijStartáµ£ = B_n*(face.neighbour-1)

       # Ïâ = cells[face.owner].var[ð.Ï]
       # Ïáµ£ = cells[face.neighbour].var[ð.Ï]

        Ïâ = face.varâ[ð.Ï]
        Ïáµ£ = face.varáµ£[ð.Ï]

        pO = cells[face.owner].var[ð.p]
        pN = cells[face.neighbour].var[ð.p]
        pâ = face.varâ[ð.p]
        páµ£ = face.varáµ£[ð.p]

       # uâ = cells[face.owner].var[ð.u]
       # uáµ£ = cells[face.neighbour].var[ð.u]
       # vâ = cells[face.owner].var[ð.v]
       # váµ£ = cells[face.neighbour].var[ð.v]

        uâ = face.varâ[ð.u]
        uáµ£ = face.varáµ£[ð.u]
        vâ = face.varâ[ð.v]
        váµ£ = face.varáµ£[ð.v]

        wâ = 0.0#cells[face.owner].var[ð.w]
        wáµ£ = 0.0#cells[face.neighbour].var[ð.w]
       #= Hââ = cells[face.owner].var[ð.Hâ]
        Hâáµ£ = cells[face.neighbour].var[ð.Hâ]
        #Î¼â = cells[face.owner].var[ð.Î¼]
        #Î¼áµ£ = cells[face.neighbour].var[ð.Î¼]
        âÏâpâ = cells[face.owner].var[ð.âÏâp]
        âÏâpáµ£ = cells[face.neighbour].var[ð.âÏâp]
        âÏâTâ = cells[face.owner].var[ð.âÏâT]
        âÏâTáµ£ = cells[face.neighbour].var[ð.âÏâT]
        âHââpâ = cells[face.owner].var[ð.âHââp]
        âHââpáµ£ = cells[face.neighbour].var[ð.âHââp]
        âHââTâ = cells[face.owner].var[ð.âHââT]
        âHââTáµ£ = cells[face.neighbour].var[ð.âHââT]
        Yââ = cells[face.owner].var[ð.Yâ]
        Yâáµ£ = cells[face.neighbour].var[ð.Yâ]
        âÏâYââ = cells[face.owner].var[ð.âÏâYâ]
        âÏâYâáµ£ = cells[face.neighbour].var[ð.âÏâYâ]
        âHââYââ = cells[face.owner].var[ð.âHââYâ]
        âHââYâáµ£ = cells[face.neighbour].var[ð.âHââYâ]
        =#
        
        Hââ = face.varâ[ð.Hâ]
        Hâáµ£ = face.varáµ£[ð.Hâ]
        âÏâpâ = face.varâ[ð.âÏâp]
        âÏâpáµ£ = face.varáµ£[ð.âÏâp]
        âÏâTâ = face.varâ[ð.âÏâT]
        âÏâTáµ£ = face.varáµ£[ð.âÏâT]
        âHââpâ = face.varâ[ð.âHââp]
        âHââpáµ£ = face.varáµ£[ð.âHââp]
        âHââTâ = face.varâ[ð.âHââT]
        âHââTáµ£ = face.varáµ£[ð.âHââT]
        Yââ = face.varâ[ð.Yâ]
        Yâáµ£ = face.varáµ£[ð.Yâ]
        âÏâYââ = face.varâ[ð.âÏâYâ]
        âÏâYâáµ£ = face.varáµ£[ð.âÏâYâ]
        âHââYââ = face.varâ[ð.âHââYâ]
        âHââYâáµ£ = face.varáµ£[ð.âHââYâ]
        câ = face.varâ[ð.c]
        cáµ£ = face.varáµ£[ð.c]

        Uââ = uâ * face.nÌ[1] + vâ * face.nÌ[2]
        Uâáµ£ = uáµ£ * face.nÌ[1] + váµ£ * face.nÌ[2]
        Uâ = 0.5 * (Uââ + Uâáµ£)
        ÎS = face.ÎS

        centerâ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centeráµ£ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ÎLR = norm(centeráµ£ - centerâ)

        ÏË¢ = 1.0 / (0.5/Ïâ + 0.5/Ïáµ£)
        #d = 0.5 * (1.0 / (Ap[face.owner]) + 1.0 / (Ap[face.neighbour]) )
        dÌ = ð.Ît / ÏË¢
        #dÌ = d / (2.0 + ÏË¢ / ð.Ît * d)
        #if d>1.e9
        #    dÌ = ð.Ît / ÏË¢
        #end
        
        # Rhie-Chow
        Uâ_RC = 0.0
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 1] * face.nÌ[1]
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 2] * face.nÌ[2]
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïâ * âÎpâx0[face.owner, 3] * face.nÌ[3]
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 1] * face.nÌ[1]
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 2] * face.nÌ[2]
        Uâ_RC += dÌ * ÏË¢ * 0.5 / Ïáµ£ * âÎpâx0[face.neighbour, 3] * face.nÌ[3]
        Uâ_RC -= dÌ * (pN-pO) / ÎLR

        RCdiffÎp = dÌ / ÎLR

        # before step
        Ïââ¿ = cells[face.owner].var[ð.Ïâ¿]
        Ïáµ£â¿ = cells[face.neighbour].var[ð.Ïâ¿]
        ÏË¢â¿ = 1.0 / (0.5/Ïââ¿ + 0.5/Ïáµ£â¿)
        Uâââ¿ = cells[face.owner].var[ð.uâ¿] * face.nÌ[1] + cells[face.owner].var[ð.vâ¿] * face.nÌ[2]
        Uâáµ£â¿ = cells[face.neighbour].var[ð.uâ¿] * face.nÌ[1] + cells[face.neighbour].var[ð.vâ¿] * face.nÌ[2]
        #Uâ_RC += dÌ * ÏË¢â¿ / ð.Ît * ( face.Uââ¿ - 0.5 * (Uâââ¿ + Uâáµ£â¿) )
        #Uâ += ( face.Uââ¿ - 0.5 * (Uâââ¿ + Uâáµ£â¿) )

        # YYL riemann
        cÌ = 0.5*(câ + cáµ£)
        Mâ = Uââ/cÌ
        Máµ£ = Uâáµ£/cÌ
        Mââº = M_func(Mâ,1.0,0.125)
        pâº = pre_func(Mâ,1.0,0.1875)
        Máµ£â» = M_func(Máµ£,-1.0,0.125)
        pâ» = pre_func(Máµ£,-1.0,0.1875)
        KLR = sqrt(0.5*(uâ^2+vâ^2+wâ^2+uáµ£^2+váµ£^2+wáµ£^2))
        Mdash = min(1.0,KLR/cÌ)
        g = -max(min(Mâ,0.0),-1.0)*min(max(Máµ£,0.0),1.0)
        Vn = (Ïâ*abs(Uââ)+Ïáµ£*abs(Uâáµ£)) / (Ïâ+Ïáµ£)
        Vnp = (1.0-g)*Vn + g*abs(Uââ)
        Vnm = (1.0-g)*Vn + g*abs(Uâáµ£)
        mÌ = 0.5*(Ïâ*(Uââ+Vnp)+Ïáµ£*(Uâáµ£-Vnm)-(1.0-Mdash)^2/cÌ*(páµ£-pâ))
        #mÌ = 0.5*(Ïâ*Uââ+Ïáµ£*Uâáµ£) + 0.5*(Ïâ*Vnp-Ïáµ£*Vnm) - 0.5*(1.0-Mdash)^2/cÌ*(páµ£-pâ)
        mÌâ = 0.5*(mÌ+abs(mÌ))
        mÌáµ£ = 0.5*(mÌ-abs(mÌ))
    
        #pâáµ£ = 0.5*(pâ+páµ£) - 0.5*KLR/cÌ*0.5*pâº*pâ»*0.5*(pâ+páµ£)/cÌ*(Uâáµ£-Uââ) + 
        #0.5*(KLR/cÌ)*0.5*(pâ+páµ£)*(pâº+pâ»-1.0) - 0.5*(pâº-pâ»)*(páµ£-pâ)

        Wâ = 0.5 * (1.0 + sign(mÌ))
        Wáµ£ = 1.0 - Wâ

        Uâ = (Wâ/Ïâ + Wáµ£/Ïáµ£)*mÌ + Uâ_RC

        WUâ = (Wâ/Ïâ + Wáµ£/Ïáµ£)*0.5*Ïâ
        WUáµ£ = (Wâ/Ïâ + Wáµ£/Ïáµ£)*0.5*Ïáµ£

        Wpâ = pâº
        Wpáµ£ = pâ»

        diffÎp = RCdiffÎp + (Wâ/Ïâ + Wáµ£/Ïáµ£)*0.5*(1.0-Mdash)^2/cÌ

        #--------------------
        # SAVE
        face.Uâ = Uâ
        #--------------------

        

        Ïâ = Wâ * Ïâ + Wáµ£ * Ïáµ£
        uâ = Wâ * uâ + Wáµ£ * uáµ£
        vâ = Wâ * vâ + Wáµ£ * váµ£
        wâ = 0.0#Wâ * wâ + Wáµ£ * wáµ£
        Hââ = Wâ * Hââ + Wáµ£ * Hâáµ£
        Yââ = Wâ * Yââ + Wáµ£ * Yâáµ£

        pâ = Wpâ*pâ + Wpáµ£*páµ£

        
        iâ = A_n*(face.owner-1)
        iáµ£ = A_n*(face.neighbour-1)


        #------------------------
        # continuity
        # p'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâpâ * Uâ * ÎS + Ïâ * diffÎp * ÎS )
        push!(A_rows, ijStartâ + 1); push!(A_cols, ijStartáµ£ + 1)
        push!(A_vals, ( Wáµ£ * âÏâpáµ£ * Uâ * ÎS - Ïâ * diffÎp * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâpáµ£ * Uâ * ÎS - Ïâ * diffÎp * ÎS )
        push!(A_rows, ijStartáµ£ + 1); push!(A_cols, ijStartâ + 1)
        push!(A_vals, -( Wâ * âÏâpâ * Uâ * ÎS + Ïâ * diffÎp * ÎS ))
        
        # u'
        iâ += 1; iáµ£ += 1
        
        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[1] * ÎS )
        push!(A_rows, ijStartâ + 1); push!(A_cols, ijStartáµ£ + 2)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[1] * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[1] * ÎS )
        push!(A_rows, ijStartáµ£ + 1); push!(A_cols, ijStartâ + 2)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[1] * ÎS ))
        
        # v'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[2] * ÎS )
        push!(A_rows, ijStartâ + 1); push!(A_cols, ijStartáµ£ + 3)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[2] * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[2] * ÎS )
        push!(A_rows, ijStartáµ£ + 1); push!(A_cols, ijStartâ + 3)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[2] * ÎS ))
        
        # T'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâTâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 1); push!(A_cols, ijStartáµ£ + 4)
        push!(A_vals, ( Wáµ£ * âÏâTáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâTáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 1); push!(A_cols, ijStartâ + 4)
        push!(A_vals, -( Wâ * âÏâTâ * Uâ * ÎS ))
        
        # Yâ'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâYââ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 1); push!(A_cols, ijStartáµ£ + 5)
        push!(A_vals, ( Wáµ£ * âÏâYâáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâYâáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 1); push!(A_cols, ijStartâ + 5)
        push!(A_vals, -( Wâ * âÏâYââ * Uâ * ÎS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâpâ * uâ * Uâ * ÎS + Wpâ * face.nÌ[1] * ÎS + Ïâ * diffÎp * uâ * ÎS )
        push!(A_rows, ijStartâ + 2); push!(A_cols, ijStartáµ£ + 1)
        push!(A_vals, ( Wáµ£ * âÏâpáµ£ * uâ * Uâ * ÎS + Wpáµ£ * face.nÌ[1] * ÎS - Ïâ * diffÎp * uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâpáµ£ * uâ * Uâ * ÎS + Wpáµ£ * face.nÌ[1] * ÎS - Ïâ * diffÎp * uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 2); push!(A_cols, ijStartâ + 1)
        push!(A_vals, -( Wâ * âÏâpâ * uâ * Uâ * ÎS + Wpâ * face.nÌ[1] * ÎS + Ïâ * diffÎp * uâ * ÎS ))
        
        # u'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[1] * uâ * ÎS + Wâ * Ïâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 2); push!(A_cols, ijStartáµ£ + 2)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[1] * uâ * ÎS + Wáµ£ * Ïáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[1] * uâ * ÎS + Wáµ£ * Ïáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 2); push!(A_cols, ijStartâ + 2)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[1] * uâ * ÎS + Wâ * Ïâ * Uâ * ÎS ))
        
        # v'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[2] * uâ * ÎS )
        push!(A_rows, ijStartâ + 2); push!(A_cols, ijStartáµ£ + 3)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[2] * uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[2] * uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 2); push!(A_cols, ijStartâ + 3)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[2] * uâ * ÎS ))

        # T'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâTâ * uâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 2); push!(A_cols, ijStartáµ£ + 4)
        push!(A_vals, ( Wáµ£ * âÏâTáµ£ * uâ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâTáµ£ * uâ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 2); push!(A_cols, ijStartâ + 4)
        push!(A_vals, -( Wâ * âÏâTâ * uâ * Uâ * ÎS ))

        # Yâ'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâYââ * uâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 2); push!(A_cols, ijStartáµ£ + 5)
        push!(A_vals, ( Wáµ£ * âÏâYâáµ£ * uâ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâYâáµ£ * uâ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 2); push!(A_cols, ijStartâ + 5)
        push!(A_vals, -( Wâ * âÏâYââ * uâ * Uâ * ÎS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] +=  ( Wâ * âÏâpâ* vâ * Uâ * ÎS + Wpâ * face.nÌ[2] * ÎS + Ïâ * diffÎp * vâ * ÎS )
        push!(A_rows, ijStartâ + 3); push!(A_cols, ijStartáµ£ + 1)
        push!(A_vals, ( Wáµ£ * âÏâpáµ£ * vâ * Uâ * ÎS + Wpáµ£ * face.nÌ[2] * ÎS - Ïâ * diffÎp * vâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâpáµ£* vâ * Uâ * ÎS + Wpáµ£ * face.nÌ[2] * ÎS - Ïâ * diffÎp * vâ * ÎS )
        push!(A_rows, ijStartáµ£ + 3); push!(A_cols, ijStartâ + 1)
        push!(A_vals, -( Wâ * âÏâpâ * vâ * Uâ * ÎS + Wpâ * face.nÌ[2] * ÎS + Ïâ * diffÎp * vâ * ÎS ))
        
        # u'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[1] * vâ * ÎS )
        push!(A_rows, ijStartâ + 3); push!(A_cols, ijStartáµ£ + 2)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[1] * vâ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[1] * vâ * ÎS )
        push!(A_rows, ijStartáµ£ + 3); push!(A_cols, ijStartâ + 2)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[1] * vâ * ÎS ))
        
        # v'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[2] * vâ * ÎS + Wâ * Ïâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 3); push!(A_cols, ijStartáµ£ + 3)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[2] * vâ * ÎS + Wáµ£ * Ïáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[2] * vâ * ÎS + Wáµ£ * Ïáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 3); push!(A_cols, ijStartâ + 3)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[2] * vâ * ÎS + Wâ * Ïâ * Uâ * ÎS ))

        # T'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâTâ * vâ * Uâ *ÎS )
        push!(A_rows, ijStartâ + 3); push!(A_cols, ijStartáµ£ + 4)
        push!(A_vals, ( Wáµ£ * âÏâTáµ£ * vâ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâTáµ£ * vâ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 3); push!(A_cols, ijStartâ + 4)
        push!(A_vals, -(  Wâ * âÏâTâ * vâ * Uâ *ÎS ))

        # Yâ'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâYââ * vâ * Uâ *ÎS )
        push!(A_rows, ijStartâ + 3); push!(A_cols, ijStartáµ£ + 5)
        push!(A_vals, ( Wáµ£ * âÏâYâáµ£ * vâ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâYâáµ£ * vâ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 3); push!(A_cols, ijStartâ + 5)
        push!(A_vals, -(  Wâ * âÏâYââ * vâ * Uâ *ÎS ))


        

        #------------------------
        # energy
        # p'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâpâ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââpâ * Uâ * ÎS + Ïâ * diffÎp * Hââ * ÎS )
        push!(A_rows, ijStartâ + 4); push!(A_cols, ijStartáµ£ + 1)
        push!(A_vals, ( Wáµ£ * âÏâpáµ£ * Hââ * Uâ * ÎS + Ïâ * Wáµ£ * âHââpáµ£ * Uâ * ÎS - Ïâ * diffÎp * Hââ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâpáµ£ * Hââ * Uâ * ÎS + Ïâ * Wáµ£ * âHââpáµ£ * Uâ * ÎS - Ïâ * diffÎp * Hââ * ÎS )
        push!(A_rows, ijStartáµ£ + 4); push!(A_cols, ijStartâ + 1)
        push!(A_vals, -( Wâ * âÏâpâ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââpâ * Uâ * ÎS + Ïâ * diffÎp * Hââ * ÎS ))
        
        # u'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[1] * Hââ * ÎS + Ïâ * Uâ * Wâ * uâ * ÎS )
        push!(A_rows, ijStartâ + 4); push!(A_cols, ijStartáµ£ + 2)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[1] * Hââ * ÎS + Ïâ * Uâ * Wáµ£ * uáµ£ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[1] * Hââ * ÎS + Ïâ * Uâ * Wáµ£ * uáµ£ * ÎS )
        push!(A_rows, ijStartáµ£ + 4); push!(A_cols, ijStartâ + 2)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[1] * Hââ * ÎS + Ïâ * Uâ * Wâ * uâ * ÎS ))

        # v'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[2] * Hââ * ÎS + Ïâ * Uâ * Wâ * vâ * ÎS )
        push!(A_rows, ijStartâ + 4); push!(A_cols, ijStartáµ£ + 3)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[2] * Hââ * ÎS + Ïâ * Uâ * Wáµ£ * váµ£ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[2] * Hââ * ÎS + Ïâ * Uâ * Wáµ£ * váµ£ * ÎS )
        push!(A_rows, ijStartáµ£ + 4); push!(A_cols, ijStartâ + 3)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[2] * Hââ * ÎS + Ïâ * Uâ * Wâ * vâ * ÎS ))

        
        # T'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâTâ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââTâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 4); push!(A_cols, ijStartáµ£ + 4)
        push!(A_vals, ( Wáµ£ * âÏâTáµ£ * Hâáµ£ * Uâ * ÎS + Ïâ * Wáµ£ * âHââTáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâTáµ£ * Hâáµ£ * Uâ * ÎS + Ïâ * Wáµ£ * âHââTáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 4); push!(A_cols, ijStartâ + 4)
        push!(A_vals, -( Wâ * âÏâTâ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââTâ * Uâ * ÎS ))

        
        # Yâ'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâYââ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââYââ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 4); push!(A_cols, ijStartáµ£ + 5)
        push!(A_vals, ( Wáµ£ * âÏâYâáµ£ * Hâáµ£ * Uâ * ÎS + Ïâ * Wáµ£ * âHââYâáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâYâáµ£ * Hâáµ£ * Uâ * ÎS + Ïâ * Wáµ£ * âHââYâáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 4); push!(A_cols, ijStartâ + 5)
        push!(A_vals, -( Wâ * âÏâYââ * Hââ * Uâ * ÎS + Ïâ * Wâ * âHââYââ * Uâ * ÎS ))
        

        

        #------------------------
        # mass fraction
        # p'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâpâ * Yââ * Uâ * ÎS + Ïâ * diffÎp * Yââ * ÎS )
        push!(A_rows, ijStartâ + 5); push!(A_cols, ijStartáµ£ + 1)
        push!(A_vals, ( Wáµ£ * âÏâpáµ£ * Yââ * Uâ * ÎS - Ïâ * diffÎp * Yââ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâpáµ£ * Yââ * Uâ * ÎS - Ïâ * diffÎp * Yââ * ÎS )
        push!(A_rows, ijStartáµ£ + 5); push!(A_cols, ijStartâ + 1)
        push!(A_vals, -( Wâ * âÏâpâ * Yââ * Uâ * ÎS + Ïâ * diffÎp * Yââ * ÎS ))
        
        # u'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[1] * Yââ * ÎS )
        push!(A_rows, ijStartâ + 5); push!(A_cols, ijStartáµ£ + 2)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[1] * Yââ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[1] * Yââ * ÎS )
        push!(A_rows, ijStartáµ£ + 5); push!(A_cols, ijStartâ + 2)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[1] * Yââ * ÎS ))

        # v'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Ïâ * WUâ * face.nÌ[2] * Yââ * ÎS )
        push!(A_rows, ijStartâ + 5); push!(A_cols, ijStartáµ£ + 3)
        push!(A_vals, ( Ïâ * WUáµ£ * face.nÌ[2] * Yââ * ÎS ))
        
        A_vals[iáµ£] -= ( Ïâ * WUáµ£ * face.nÌ[2] * Yââ * ÎS )
        push!(A_rows, ijStartáµ£ + 5); push!(A_cols, ijStartâ + 3)
        push!(A_vals, -( Ïâ * WUâ * face.nÌ[2] * Yââ * ÎS ))

        
        # T'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâTâ * Yââ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 5); push!(A_cols, ijStartáµ£ + 4)
        push!(A_vals, ( Wáµ£ * âÏâTáµ£ * Yââ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâTáµ£ * Yââ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 5); push!(A_cols, ijStartâ + 4)
        push!(A_vals, -( Wâ * âÏâTâ * Yââ * Uâ * ÎS ))

        
        # Yâ'
        iâ += 1; iáµ£ += 1

        A_vals[iâ] += ( Wâ * âÏâYââ * Yââ * Uâ * ÎS + Ïâ * Wâ * Uâ * ÎS )
        push!(A_rows, ijStartâ + 5); push!(A_cols, ijStartáµ£ + 5)
        push!(A_vals, ( Wáµ£ * âÏâYâáµ£ * Yââ * Uâ * ÎS + Ïâ * Wáµ£ * Uâ * ÎS ))
        
        A_vals[iáµ£] -= ( Wáµ£ * âÏâYâáµ£ * Yââ * Uâ * ÎS + Ïâ * Wáµ£ * Uâ * ÎS )
        push!(A_rows, ijStartáµ£ + 5); push!(A_cols, ijStartâ + 5)
        push!(A_vals, -( Wâ * âÏâYââ * Yââ * Uâ * ÎS + Ïâ * Wâ * Uâ * ÎS ))
        

        # ----------------------------

        B[ijStartâ + 1] -= ( Ïâ * Uâ * ÎS )
        B[ijStartáµ£ + 1] += ( Ïâ * Uâ * ÎS )

        B[ijStartâ + 2] -= ( Ïâ * uâ * Uâ * ÎS + pâ * face.nÌ[1] * ÎS )
        B[ijStartáµ£ + 2] += ( Ïâ * uâ * Uâ * ÎS + pâ * face.nÌ[1] * ÎS )

        B[ijStartâ + 3] -= ( Ïâ * vâ * Uâ * ÎS + pâ * face.nÌ[2] * ÎS )
        B[ijStartáµ£ + 3] += ( Ïâ * vâ * Uâ * ÎS + pâ * face.nÌ[2] * ÎS )

        B[ijStartâ + 4] -= ( Ïâ * Hââ * Uâ * ÎS )
        B[ijStartáµ£ + 4] += ( Ïâ * Hââ * Uâ * ÎS )

        B[ijStartâ + 5] -= ( Ïâ * Yââ * Uâ * ÎS )
        B[ijStartáµ£ + 5] += ( Ïâ * Yââ * Uâ * ÎS )




    end


    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    #bc_wall = []
    #append!( bc_wall, faces_boundary_top )
    #append!( bc_wall, faces_boundary_bottom )
    #append!( bc_wall, faces_boundary_left )
    #append!( bc_wall, faces_boundary_right )


    coupled_boundary!(
        ð,cells,faces,
        faces_boundary_top, 
        ð.top_p_BCtype, ð.top_p_BCValue, 
        ð.top_u_BCtype, ð.top_u_BCValue, 
        ð.top_v_BCtype, ð.top_v_BCValue, 
        ð.top_T_BCtype, ð.top_T_BCValue, 
        ð.top_Y_BCtype, ð.top_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(ð,cells,faces,
        faces_boundary_bottom, 
        ð.bottom_p_BCtype, ð.bottom_p_BCValue, 
        ð.bottom_u_BCtype, ð.bottom_u_BCValue, 
        ð.bottom_v_BCtype, ð.bottom_v_BCValue, 
        ð.bottom_T_BCtype, ð.bottom_T_BCValue, 
        ð.bottom_Y_BCtype, ð.bottom_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(ð,cells,faces,
        faces_boundary_left, 
        ð.left_p_BCtype, ð.left_p_BCValue, 
        ð.left_u_BCtype, ð.left_u_BCValue, 
        ð.left_v_BCtype, ð.left_v_BCValue, 
        ð.left_T_BCtype, ð.left_T_BCValue, 
        ð.left_Y_BCtype, ð.left_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(ð,cells,faces,
        faces_boundary_right, 
        ð.right_p_BCtype, ð.right_p_BCValue, 
        ð.right_u_BCtype, ð.right_u_BCValue, 
        ð.right_v_BCtype, ð.right_v_BCValue, 
        ð.right_T_BCtype, ð.right_T_BCValue, 
        ð.right_Y_BCtype, ð.right_Y_BCValue,
        B_n, A_n, A_vals, B)


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ÎQ = solve(ps, A, B)

    #ÎQ = A\B

    #ml = ruge_stuben(A)
    #ÎQ = solve(ml, A)
    #P = aspreconditioner(ml)
    #ÎQ = bicgstabl(A, B, Pl = P)
    #ÎQ = gmres(A, B)

    #ÎQ = A\B

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

        cell.var[ð.p] += relax_p * ÎQ[ijStart + 1]
        cell.var[ð.u] += relax_U * ÎQ[ijStart + 2]
        cell.var[ð.v] += relax_U * ÎQ[ijStart + 3]
        cell.var[ð.T] += relax_T * ÎQ[ijStart + 4]
        cell.var[ð.Yâ] += relax_Y * ÎQ[ijStart + 5]

        cell.var[ð.p] = max(cell.var[ð.p],10.0)
        cell.var[ð.T] = max(cell.var[ð.T],10.0)
        
        norm_total += cell.var[ð.p]^2
        norm_total += cell.var[ð.u]^2
        norm_total += cell.var[ð.v]^2
        norm_total += cell.var[ð.T]^2
        norm_total += cell.var[ð.Yâ]^2

        diagon += 1
    end

    norm_total = sqrt(norm_total)

    #sleep(1000.0)


    return norm(ÎQ)/length(cells)/5
   


end



function coupled_boundary!(
    ð::controls,
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
        
        ijStartâ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        Ïâ = cells[face.owner].var[ð.Ï]
        âÏâpâ = cells[face.owner].var[ð.âÏâp]
        âÏâTâ = cells[face.owner].var[ð.âÏâT]
        âHââpâ = cells[face.owner].var[ð.âHââp]
        âHââTâ = cells[face.owner].var[ð.âHââT]
        âHââYââ = cells[face.owner].var[ð.âHââYâ]
        âÏâYââ = cells[face.owner].var[ð.âÏâYâ]
        pâ = cells[face.owner].var[ð.p]
        Hââ = cells[face.owner].var[ð.Hâ]
        Yââ = cells[face.owner].var[ð.Yâ]
        Tâ = cells[face.owner].var[ð.T]

        ÎS = face.ÎS

        Uâ = 0.0
        Uâ += cells[face.owner].var[ð.u]*face.nÌ[1]
        Uâ += cells[face.owner].var[ð.v]*face.nÌ[2]
        Uâ += cells[face.owner].var[ð.w]*face.nÌ[3]

        uâ = cells[face.owner].var[ð.u] - Uâ * face.nÌ[1]
        vâ = cells[face.owner].var[ð.v] - Uâ * face.nÌ[2]
        wâ = 0.0#cells[face.owner].var[ð.w] - Uâ * face.nÌ[3]

        Uâ = uâ * face.nÌ[1] + vâ * face.nÌ[2] + wâ * face.nÌ[3]

        id = []
        push!(id,i)
        push!(id,i+5)
        push!(id,i+10)
        push!(id,i+15)
        push!(id,i+20)

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pâ = cells[face.owner].var[ð.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pâ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pâ = p_BCValue(ð.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uâ = cells[face.owner].var[ð.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uâ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uâ = uâ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uâ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uâ = u_BCValue(ð.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vâ = cells[face.owner].var[ð.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vâ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vâ = vâ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vâ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vâ = v_BCValue(ð.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tâ = cells[face.owner].var[ð.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tâ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tâ = T_BCValue(ð.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yââ = cells[face.owner].var[ð.Yâ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yââ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yââ = Y_BCValue(ð.time)
        end
        
        Uâ = uâ * face.nÌ[1] + vâ * face.nÌ[2] + wâ * face.nÌ[3]

        Ïâ, Hââ, câ = faceEOS!(ð,pâ,uâ,vâ,wâ,Tâ,Yââ)
        
        # continuity
        i += 1
        A_vals[i] += coeff_p * (âÏâpâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_u * (Ïâ * face.nÌ[1] * ÎS)
        i += 1
        A_vals[i] += coeff_v * (Ïâ * face.nÌ[2] * ÎS)
        i += 1
        A_vals[i] += coeff_T * (âÏâTâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_Y * (âÏâYââ * Uâ * ÎS)

        
        # x-momentum
        i += 1
        A_vals[i] += coeff_p * (âÏâpâ * uâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_u * (Ïâ * Uâ * ÎS + Ïâ * uâ * face.nÌ[1] * ÎS)
        i += 1
        A_vals[i] += coeff_v * (Ïâ * uâ * face.nÌ[2] * ÎS)
        i += 1
        A_vals[i] += coeff_T * (âÏâTâ * uâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_Y * (âÏâYââ * uâ * Uâ * ÎS)

        
        # y-momentum
        i += 1
        A_vals[i] += coeff_p * (âÏâpâ * vâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_u * (Ïâ * vâ * face.nÌ[1] * ÎS)
        i += 1
        A_vals[i] += coeff_v * (Ïâ * Uâ * ÎS + Ïâ * vâ * face.nÌ[2] * ÎS)
        i += 1
        A_vals[i] += coeff_T * (âÏâTâ * vâ * Uâ * ÎS)
        i += 1
        A_vals[i] += coeff_Y * (âÏâYââ * vâ * Uâ * ÎS)


        # energy
        i += 1
        A_vals[i] += coeff_p * (âÏâpâ * Uâ * Hââ * ÎS + Ïâ * Uâ * âHââpâ * ÎS)
        i += 1
        A_vals[i] += coeff_u * (Ïâ * face.nÌ[1] * Hââ * ÎS + Ïâ * Uâ * uâ * ÎS)
        i += 1
        A_vals[i] += coeff_v * (Ïâ * face.nÌ[2] * Hââ * ÎS + Ïâ * Uâ * vâ * ÎS)
        i += 1
        A_vals[i] += coeff_T * (âÏâTâ * Uâ * Hââ * ÎS + Ïâ * Uâ * âHââTâ * ÎS)
        i += 1
        A_vals[i] += coeff_Y * (âÏâYââ * Uâ * Hââ * ÎS + Ïâ * Uâ * âHââYââ * ÎS)


        # massfraction
        i += 1
        A_vals[i] += coeff_p * (âÏâpâ * Uâ * Yââ * ÎS)# + Ïâ * Hââ * ð.Ît/Ïâ / ÎLR * ÎS
        i += 1
        A_vals[i] += coeff_u * (Ïâ * face.nÌ[1] * Yââ * ÎS)
        i += 1
        A_vals[i] += coeff_v * (Ïâ * face.nÌ[2] * Yââ * ÎS)
        i += 1
        A_vals[i] += coeff_T * (âÏâTâ * Uâ * Yââ * ÎS)
        i += 1
        A_vals[i] += coeff_Y * (âÏâYââ * Uâ * Yââ * ÎS + Ïâ * Uâ * ÎS)


        B[ijStartâ + 1] -= ( Ïâ * Uâ * ÎS )
        B[ijStartâ + 2] -= ( Ïâ * uâ * Uâ * ÎS + pâ * face.nÌ[1] * ÎS )
        B[ijStartâ + 3] -= ( Ïâ * vâ * Uâ * ÎS + pâ * face.nÌ[2] * ÎS )
        B[ijStartâ + 4] -= ( Ïâ * Hââ * Uâ * ÎS )
        B[ijStartâ + 5] -= ( Ïâ * Yââ * Uâ * ÎS )
        

    end
 


end












function coupled_Ap_boundary!(
    ð::controls,
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

        Ïâ = cells[face.owner].var[ð.Ï]
        pâ = cells[face.owner].var[ð.p]
        Hââ = cells[face.owner].var[ð.Hâ]
        Yââ = cells[face.owner].var[ð.Yâ]
        Tâ = cells[face.owner].var[ð.T]

        ÎS = face.ÎS

        Uâ = 0.0
        Uâ += cells[face.owner].var[ð.u]*face.nÌ[1]
        Uâ += cells[face.owner].var[ð.v]*face.nÌ[2]
        Uâ += cells[face.owner].var[ð.w]*face.nÌ[3]

        uâ = cells[face.owner].var[ð.u] - Uâ * face.nÌ[1]
        vâ = cells[face.owner].var[ð.v] - Uâ * face.nÌ[2]
        wâ = 0.0#cells[face.owner].var[ð.w] - Uâ * face.nÌ[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pâ = cells[face.owner].var[ð.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pâ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pâ = p_BCValue(ð.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uâ = cells[face.owner].var[ð.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uâ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uâ = uâ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uâ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uâ = u_BCValue(ð.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vâ = cells[face.owner].var[ð.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vâ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vâ = vâ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vâ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vâ = v_BCValue(ð.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tâ = cells[face.owner].var[ð.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tâ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tâ = T_BCValue(ð.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yââ = cells[face.owner].var[ð.Yâ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yââ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yââ = Y_BCValue(ð.time)
        end
        
        Uâ = uâ * face.nÌ[1] + vâ * face.nÌ[2] + wâ * face.nÌ[3]

        Ïâ, Hââ, câ = faceEOS!(ð,pâ,uâ,vâ,wâ,Tâ,Yââ)
        

        flux = Ïâ * Uâ * ÎS
        Ap[face.owner] += flux / cells[face.owner].Î©


    end
 


end


function M_func(M::Float64, op::Float64, Î±::Float64)
    mu = 0.0
	if abs(M) > 1.0 
		mu = 0.5*(M + op*abs(M))
	else
		mu = op*0.25*(M + op)^2.0 + op*Î±*(M*M-1.0)^2.0
    end
	
	return mu
end

function pre_func(M::Float64, op::Float64, Î±::Float64)
    mu = 0.0
	if abs(M) > 1.0
		mu = 0.5*(1.0 + op*sign(M) )
	else
		mu = 0.25*(M + op)^2.0*(2.0-op*M) + op*Î±*M*(M*M-1.0)^2.0
    end
	
	return mu;
end


