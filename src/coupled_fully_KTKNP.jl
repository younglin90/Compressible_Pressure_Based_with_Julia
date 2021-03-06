

function coupled!(
    ๐::controls,
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

        ฮฉ = cell.ฮฉ
        ฮt = ๐.ฮt
        u = cell.var[๐.u]
        v = cell.var[๐.v]
        ฯ = cell.var[๐.ฯ]
        Hโ = cell.var[๐.Hโ]
        p = cell.var[๐.p]
        Yโ = cell.var[๐.Yโ]
        โHโโp = cell.var[๐.โHโโp]
        โHโโT = cell.var[๐.โHโโT]
        โฯโp = cell.var[๐.โฯโp]
        โฯโT = cell.var[๐.โฯโT]
        โฯโYโ = cell.var[๐.โฯโYโ]
        โHโโYโ = cell.var[๐.โHโโYโ]
        ฯโฟ = cell.var[๐.ฯโฟ]
        uโฟ = cell.var[๐.uโฟ]
        vโฟ = cell.var[๐.vโฟ]
        Hโโฟ = cell.var[๐.Hโโฟ]
        pโฟ = cell.var[๐.pโฟ]
        Yโโฟ = cell.var[๐.Yโโฟ]
        ฯโฟโปยน = cell.var[๐.ฯโฟโปยน]
        uโฟโปยน = cell.var[๐.uโฟโปยน]
        vโฟโปยน = cell.var[๐.vโฟโปยน]
        Hโโฟโปยน = cell.var[๐.Hโโฟโปยน]
        pโฟโปยน = cell.var[๐.pโฟโปยน]
        Yโโฟโปยน = cell.var[๐.Yโโฟโปยน]

        c_ฮt = 1.0
        if ๐.temporal_discretizationScheme == "1st"  
            c_ฮt = 1.0
        elseif ๐.temporal_discretizationScheme == "2nd"
            c_ฮt = 1.5
        end
        #println(pโฟ,uโฟ,vโฟ,Hโโฟ)
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = c_ฮt * ( โฯโp*ฮฉ/ฮt )
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = c_ฮt * ( โฯโT*ฮฉ/ฮt )
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 5
        A_vals[i] = c_ฮt * ( โฯโYโ*ฮฉ/ฮt )

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = c_ฮt * ( โฯโp*ฮฉ/ฮt * u ) - โฯโp*๐.gravity[1]*ฮฉ

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = c_ฮt * ( ฯ*ฮฉ/ฮt )
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = c_ฮt * ( โฯโT*ฮฉ/ฮt * u ) - โฯโT*๐.gravity[1]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 5
        A_vals[i] = c_ฮt * ( โฯโYโ*ฮฉ/ฮt * u ) - โฯโYโ*๐.gravity[1]*ฮฉ


        # y-momentum
        #g = -9.8
        #g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = c_ฮt * ( โฯโp*ฮฉ/ฮt * v ) - โฯโp*๐.gravity[2]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = c_ฮt * ( ฯ*ฮฉ/ฮt )

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = c_ฮt * ( โฯโT*ฮฉ/ฮt * v ) - โฯโT*๐.gravity[2]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 5
        A_vals[i] = c_ฮt * ( โฯโYโ*ฮฉ/ฮt * v ) - โฯโYโ*๐.gravity[2]*ฮฉ




        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = c_ฮt * ( โฯโp*ฮฉ/ฮt * Hโ + โHโโp*ฮฉ/ฮt * ฯ - ฮฉ/ฮt )

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = c_ฮt * ( u*ฮฉ/ฮt * ฯ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = c_ฮt * ( v*ฮฉ/ฮt * ฯ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = c_ฮt * ( โฯโT*ฮฉ/ฮt * Hโ + โHโโT*ฮฉ/ฮt * ฯ )
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 5
        A_vals[i] = c_ฮt * ( โฯโYโ*ฮฉ/ฮt * Hโ + โHโโYโ*ฮฉ/ฮt * ฯ )



        # mass fraction
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 1
        A_vals[i] = c_ฮt * ( โฯโp*ฮฉ/ฮt * Yโ )

        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 4
        A_vals[i] = c_ฮt * ( โฯโT*ฮฉ/ฮt * Yโ )
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 5
        A_vals[i] = c_ฮt * ( โฯโYโ*ฮฉ/ฮt * Yโ + ฮฉ/ฮt * ฯ )
        
        
        
        # B
        if ๐.temporal_discretizationScheme == "1st"
            B[ijStart + 1] = -(ฯ - ฯโฟ)*ฮฉ/ฮt
            B[ijStart + 2] = -(ฯ*u - ฯโฟ*uโฟ)*ฮฉ/ฮt + ฯ*๐.gravity[1]*ฮฉ 
            B[ijStart + 3] = -(ฯ*v - ฯโฟ*vโฟ)*ฮฉ/ฮt + ฯ*๐.gravity[2]*ฮฉ 
            B[ijStart + 4] = -(ฯ*Hโ - ฯโฟ*Hโโฟ)*ฮฉ/ฮt + (p - pโฟ)*ฮฉ/ฮt
            B[ijStart + 5] = -(ฯ*Yโ - ฯโฟ*Yโโฟ)*ฮฉ/ฮt
        elseif ๐.temporal_discretizationScheme == "2nd"
            B[ijStart + 1] = -(1.5*ฯ - 2.0*ฯโฟ + 0.5*ฯโฟโปยน)*ฮฉ/ฮt
            B[ijStart + 2] = -(1.5*ฯ*u - 2.0*ฯโฟ*uโฟ + 0.5*ฯโฟโปยน*uโฟโปยน)*ฮฉ/ฮt + ฯ*๐.gravity[1]*ฮฉ 
            B[ijStart + 3] = -(1.5*ฯ*v - 2.0*ฯโฟ*vโฟ + 0.5*ฯโฟโปยน*vโฟโปยน)*ฮฉ/ฮt + ฯ*๐.gravity[2]*ฮฉ 
            B[ijStart + 4] = -(1.5*ฯ*Hโ - 2.0*ฯโฟ*Hโโฟ + 0.5*ฯโฟโปยน*Hโโฟโปยน)*ฮฉ/ฮt + (1.5*p - 2.0*pโฟ + 0.5*pโฟโปยน)*ฮฉ/ฮt
            B[ijStart + 5] = -(1.5*ฯ*Yโ - 2.0*ฯโฟ*Yโโฟ + 0.5*ฯโฟโปยน*Yโโฟโปยน)*ฮฉ/ฮt
        end


        diagon += 1


    end

    




    โฮpโx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pโ = 0.5 * (cells[face.owner].var[๐.p] + cells[face.neighbour].var[๐.p])
        โฮpโx0[face.owner, 1] += pโ * face.nฬ[1] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 2] += pโ * face.nฬ[2] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 3] += pโ * face.nฬ[3] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.neighbour, 1] -= pโ * face.nฬ[1] * face.ฮS / cells[face.neighbour].ฮฉ
        โฮpโx0[face.neighbour, 2] -= pโ * face.nฬ[2] * face.ฮS / cells[face.neighbour].ฮฉ
        โฮpโx0[face.neighbour, 3] -= pโ * face.nฬ[3] * face.ฮS / cells[face.neighbour].ฮฉ
    end

    for face in faces_boundary
        pโ = cells[face.owner].var[๐.p]
        โฮpโx0[face.owner, 1] += pโ * face.nฬ[1] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 2] += pโ * face.nฬ[2] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 3] += pโ * face.nฬ[3] * face.ฮS / cells[face.owner].ฮฉ
    end




    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartโ = B_n*(face.owner-1)
        ijStartแตฃ = B_n*(face.neighbour-1)

       # ฯโ = cells[face.owner].var[๐.ฯ]
       # ฯแตฃ = cells[face.neighbour].var[๐.ฯ]

        ฯโ = face.varโ[๐.ฯ]
        ฯแตฃ = face.varแตฃ[๐.ฯ]

        pO = cells[face.owner].var[๐.p]
        pN = cells[face.neighbour].var[๐.p]
        pโ = face.varโ[๐.p]
        pแตฃ = face.varแตฃ[๐.p]

       # uโ = cells[face.owner].var[๐.u]
       # uแตฃ = cells[face.neighbour].var[๐.u]
       # vโ = cells[face.owner].var[๐.v]
       # vแตฃ = cells[face.neighbour].var[๐.v]

        uโ = face.varโ[๐.u]
        uแตฃ = face.varแตฃ[๐.u]
        vโ = face.varโ[๐.v]
        vแตฃ = face.varแตฃ[๐.v]

        wโ = 0.0#cells[face.owner].var[๐.w]
        wแตฃ = 0.0#cells[face.neighbour].var[๐.w]
       #= Hโโ = cells[face.owner].var[๐.Hโ]
        Hโแตฃ = cells[face.neighbour].var[๐.Hโ]
        #ฮผโ = cells[face.owner].var[๐.ฮผ]
        #ฮผแตฃ = cells[face.neighbour].var[๐.ฮผ]
        โฯโpโ = cells[face.owner].var[๐.โฯโp]
        โฯโpแตฃ = cells[face.neighbour].var[๐.โฯโp]
        โฯโTโ = cells[face.owner].var[๐.โฯโT]
        โฯโTแตฃ = cells[face.neighbour].var[๐.โฯโT]
        โHโโpโ = cells[face.owner].var[๐.โHโโp]
        โHโโpแตฃ = cells[face.neighbour].var[๐.โHโโp]
        โHโโTโ = cells[face.owner].var[๐.โHโโT]
        โHโโTแตฃ = cells[face.neighbour].var[๐.โHโโT]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Yโแตฃ = cells[face.neighbour].var[๐.Yโ]
        โฯโYโโ = cells[face.owner].var[๐.โฯโYโ]
        โฯโYโแตฃ = cells[face.neighbour].var[๐.โฯโYโ]
        โHโโYโโ = cells[face.owner].var[๐.โHโโYโ]
        โHโโYโแตฃ = cells[face.neighbour].var[๐.โHโโYโ]
        =#
        
        Hโโ = face.varโ[๐.Hโ]
        Hโแตฃ = face.varแตฃ[๐.Hโ]
        โฯโpโ = face.varโ[๐.โฯโp]
        โฯโpแตฃ = face.varแตฃ[๐.โฯโp]
        โฯโTโ = face.varโ[๐.โฯโT]
        โฯโTแตฃ = face.varแตฃ[๐.โฯโT]
        โHโโpโ = face.varโ[๐.โHโโp]
        โHโโpแตฃ = face.varแตฃ[๐.โHโโp]
        โHโโTโ = face.varโ[๐.โHโโT]
        โHโโTแตฃ = face.varแตฃ[๐.โHโโT]
        Yโโ = face.varโ[๐.Yโ]
        Yโแตฃ = face.varแตฃ[๐.Yโ]
        โฯโYโโ = face.varโ[๐.โฯโYโ]
        โฯโYโแตฃ = face.varแตฃ[๐.โฯโYโ]
        โHโโYโโ = face.varโ[๐.โHโโYโ]
        โHโโYโแตฃ = face.varแตฃ[๐.โHโโYโ]
        cโ = face.varโ[๐.c]
        cแตฃ = face.varแตฃ[๐.c]

        Uโโ = uโ * face.nฬ[1] + vโ * face.nฬ[2]
        Uโแตฃ = uแตฃ * face.nฬ[1] + vแตฃ * face.nฬ[2]
        Uโ = 0.5 * (Uโโ + Uโแตฃ)
        ฮS = face.ฮS

        centerโ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerแตฃ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ฮLR = norm(centerแตฃ - centerโ)

        ฯหข = 1.0 / (0.5/ฯโ + 0.5/ฯแตฃ)
        #d = 0.5 * (1.0 / (Ap[face.owner]) + 1.0 / (Ap[face.neighbour]) )
        dฬ = ๐.ฮt / ฯหข
        #dฬ = d / (2.0 + ฯหข / ๐.ฮt * d)
        #if d>1.e9
        #    dฬ = ๐.ฮt / ฯหข
        #end
        
        # Rhie-Chow
        Uโ_RC = 0.0
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 1] * face.nฬ[1]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 2] * face.nฬ[2]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 3] * face.nฬ[3]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 1] * face.nฬ[1]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 2] * face.nฬ[2]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 3] * face.nฬ[3]
        Uโ_RC -= dฬ * (pN-pO) / ฮLR

        RCdiffฮp = dฬ / ฮLR

        # before step
        ฯโโฟ = cells[face.owner].var[๐.ฯโฟ]
        ฯแตฃโฟ = cells[face.neighbour].var[๐.ฯโฟ]
        ฯหขโฟ = 1.0 / (0.5/ฯโโฟ + 0.5/ฯแตฃโฟ)
        Uโโโฟ = cells[face.owner].var[๐.uโฟ] * face.nฬ[1] + cells[face.owner].var[๐.vโฟ] * face.nฬ[2]
        Uโแตฃโฟ = cells[face.neighbour].var[๐.uโฟ] * face.nฬ[1] + cells[face.neighbour].var[๐.vโฟ] * face.nฬ[2]
        #Uโ_RC += dฬ * ฯหขโฟ / ๐.ฮt * ( face.Uโโฟ - 0.5 * (Uโโโฟ + Uโแตฃโฟ) )
        #Uโ += ( face.Uโโฟ - 0.5 * (Uโโโฟ + Uโแตฃโฟ) )

        # YYL riemann
        cฬ = 0.5*(cโ + cแตฃ)
        Mโ = Uโโ/cฬ
        Mแตฃ = Uโแตฃ/cฬ
        # calculate M+ and P+ for left state
        Mโโบ = M_func(Mโ,1.0,0.125)
        pโบ = pre_func(Mโ,1.0,0.1875)
        # calculate M- and P- for left state
        Mแตฃโป = M_func(Mแตฃ,-1.0,0.125)
        pโป = pre_func(Mแตฃ,-1.0,0.1875)
        KLR = sqrt(0.5*(uโ^2+vโ^2+wโ^2+uแตฃ^2+vแตฃ^2+wแตฃ^2))
        Mdash = min(1.0,KLR/cฬ)

        Cdiffฮp = 0.5*(1.0-Mdash)^2/cฬ /(0.5*(ฯโ+ฯแตฃ))
        WUโ = Mโโบ/Mโ
        if abs(Mโ) <= 0.00001
            WUโ = 0.5
        end
        WUแตฃ = Mแตฃโป/Mแตฃ
        if abs(Mแตฃ) <= 0.00001
            WUแตฃ = 0.5
        end

	    #Uโ = WUโ*Uโโ + WUแตฃ*Uโแตฃ + Uโ_RC - Cdiffฮp * (pแตฃ-pโ)
	    #Uโ = WUโ*Uโโ + WUแตฃ*Uโแตฃ - Cdiffฮp * (pแตฃ-pโ)


        Fmax =  max(max(Uโโ+cโ,Uโแตฃ+cแตฃ),0.0)
        Fmin = -min(min(Uโโ-cโ,Uโแตฃ-cแตฃ),0.0)
    
        ฮฑโบ = Fmax / (Fmax+Fmin)
        ฮฑโป = Fmin / (Fmax+Fmin)
        ฮฑโบโป = Fmax*Fmin / (Fmax+Fmin)
    
    


        Wpโ = pโบ
        Wpแตฃ = pโป

        #--------------------
        # SAVE
        face.Uโ = Uโ
        #--------------------


        ฯโ = (ฮฑโบ * Uโโ + ฮฑโบโป)
        ฯแตฃ = (ฮฑโป * Uโแตฃ - ฮฑโบโป)
        ฯโแตฃ = ฯโ+ฯแตฃ
        
        #Wโ = 1.0 #0.5 * (1.0 + sign(Uโ))
        #Wแตฃ = 1.0 #1.0 - Wโ

        #KLRโ = sqrt(uโ^2+vโ^2+wโ^2)
        #KLRแตฃ = sqrt(uแตฃ^2+vแตฃ^2+wแตฃ^2)
        #KLR = max(KLRโ,KLRแตฃ)
        #cmax = min(cโ,cแตฃ)
        #Mdash = min(1.0,KLR/cmax)
        ฮf = 1.0#min(1.0, 500.0*Mdash)
        #ฮf = sin(0.5*pi*Mdash)

        Wโ = (ฮf) * ( 1.0 ) + (1.0-ฮf) * ( 0.5 * (1.0 + sign(ฯโแตฃ)) )
        Wแตฃ = (ฮf) * ( 1.0 ) + (1.0-ฮf) * ( 0.5 * (1.0 - sign(ฯโแตฃ)) )
        
        ฯโบ = (ฮf) * ( ฯโ ) + (1.0-ฮf) * ( ฯโแตฃ + Uโ_RC ) + Uโ_RC
        ฯโป = (ฮf) * ( ฯแตฃ ) + (1.0-ฮf) * ( ฯโแตฃ + Uโ_RC ) + Uโ_RC
        
        CโUโ = (ฮf) * ( ฮฑโบ ) + (1.0-ฮf) * ( ฮฑโบ + ฮฑโป )
        CโUแตฃ = (ฮf) * ( ฮฑโป ) + (1.0-ฮf) * ( ฮฑโบ + ฮฑโป )
        
        diffฮp = (1.0-ฮf) * RCdiffฮp + RCdiffฮp
        
        ฯโบ = Wโ * ฯโ + (1.0-ฮf) * ( Wแตฃ * ฯแตฃ )
        uโบ = Wโ * uโ + (1.0-ฮf) * ( Wแตฃ * uแตฃ )
        vโบ = Wโ * vโ + (1.0-ฮf) * ( Wแตฃ * vแตฃ )
        Hโโบ = Wโ * Hโโ + (1.0-ฮf) * ( Wแตฃ * Hโแตฃ )
        Yโโบ = Wโ * Yโโ + (1.0-ฮf) * ( Wแตฃ * Yโแตฃ )
        
        ฯโป = Wแตฃ * ฯแตฃ + (1.0-ฮf) * ( Wโ * ฯโ )
        uโป = Wแตฃ * uแตฃ + (1.0-ฮf) * ( Wโ * uโ )
        vโป = Wแตฃ * vแตฃ + (1.0-ฮf) * ( Wโ * vโ )
        Hโโป = Wแตฃ * Hโแตฃ + (1.0-ฮf) * ( Wโ * Hโโ )
        Yโโป = Wแตฃ * Yโแตฃ + (1.0-ฮf) * ( Wโ * Yโโ )

        pโ = Wpโ*pโ + Wpแตฃ*pแตฃ

        
        iโ = A_n*(face.owner-1)
        iแตฃ = A_n*(face.neighbour-1)


        #------------------------
        # continuity
        # p'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโpโ * ฯโบ * ฮS + ฯโบ * diffฮp * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ * ฯโป * ฮS - ฯโป * diffฮp * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * ฯโป * ฮS - ฯโป * diffฮp * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ * ฯโบ * ฮS + ฯโบ * diffฮp * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1
        
        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[1] * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[1] * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[1] * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[1] * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[2] * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[2] * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[2] * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[2] * ฮS ))
        
        # T'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโTโ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ * ฯโบ * ฮS ))
        
        # Yโ'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโYโโ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * ฯโบ * ฮS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโpโ * uโบ * ฯโบ * ฮS + Wpโ * face.nฬ[1] * ฮS + ฯโบ * diffฮp * uโบ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ * uโป * ฯโป * ฮS + Wpแตฃ * face.nฬ[1] * ฮS - ฯโป * diffฮp * uโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * uโป * ฯโป * ฮS + Wpแตฃ * face.nฬ[1] * ฮS - ฯโป * diffฮp * uโป * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ * uโบ * ฯโบ * ฮS + Wpโ * face.nฬ[1] * ฮS + ฯโบ * diffฮp * uโบ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[1] * uโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[1] * uโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[1] * uโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[1] * uโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[2] * uโบ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[2] * uโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[2] * uโป * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[2] * uโบ * ฮS ))

        # T'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโTโ * uโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ * uโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * uโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ * uโบ * ฯโบ * ฮS ))

        # Yโ'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโYโโ * uโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * uโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * uโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * uโบ * ฯโบ * ฮS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] +=  ( Wโ * โฯโpโ * vโบ * ฯโบ * ฮS + Wpโ * face.nฬ[2] * ฮS + ฯโบ * diffฮp * vโบ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ * vโป * ฯโป * ฮS + Wpแตฃ * face.nฬ[2] * ฮS - ฯโป * diffฮp * vโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * vโป * ฯโป * ฮS + Wpแตฃ * face.nฬ[2] * ฮS - ฯโป * diffฮp * vโป * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ * vโบ * ฯโบ * ฮS + Wpโ * face.nฬ[2] * ฮS + ฯโบ * diffฮp * vโบ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[1] * vโบ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[1] * vโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[1] * vโป * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[1] * vโบ * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[2] * vโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[2] * vโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[2] * vโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[2] * vโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS ))

        # T'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโTโ * vโบ * ฯโบ *ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ * vโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * vโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ * vโบ * ฯโบ *ฮS ))

        # Yโ'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโYโโ * vโบ * ฯโบ *ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * vโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * vโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * vโบ * ฯโบ *ฮS ))


        

        #------------------------
        # energy
        # p'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโpโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโpโ * ฯโบ * ฮS + ฯโบ * diffฮp * Hโโบ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโpแตฃ * ฯโป * ฮS - ฯโป * diffฮp * Hโโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโpแตฃ * ฯโป * ฮS - ฯโป * diffฮp * Hโโป * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโpโ * ฯโบ * ฮS + ฯโบ * diffฮp * Hโโบ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[1] * Hโโบ * ฮS + Wโ * ฯโบ * ฯโบ * uโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[1] * Hโโป * ฮS + Wแตฃ * ฯโป * ฯโป * uแตฃ * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[1] * Hโโป * ฮS + Wแตฃ * ฯโป * ฯโป * uแตฃ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[1] * Hโโบ * ฮS + Wโ * ฯโบ * ฯโบ * uโ * ฮS ))

        # v'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[2] * Hโโบ * ฮS + Wโ * ฯโบ * ฯโบ * vโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[2] * Hโโป * ฮS + Wแตฃ * ฯโป * ฯโป * vแตฃ * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[2] * Hโโป * ฮS + Wแตฃ * ฯโป * ฯโป * vแตฃ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[2] * Hโโบ * ฮS + Wโ * ฯโบ * ฯโบ * vโ * ฮS ))

        
        # T'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโTโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโTโ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโTแตฃ * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโTแตฃ * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโTโ * ฯโบ * ฮS ))

        
        # Yโ'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโYโโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโYโโ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโYโแตฃ * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * Hโโป * ฯโป * ฮS + Wแตฃ * ฯโป * โHโโYโแตฃ * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * Hโโบ * ฯโบ * ฮS + Wโ * ฯโบ * โHโโYโโ * ฯโบ * ฮS ))
        
        

        #------------------------
        # mass fraction
        # p'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโpโ * Yโโบ * ฯโบ * ฮS + ฯโบ * diffฮp * Yโโบ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ * Yโโป * ฯโป * ฮS - ฯโป * diffฮp * Yโโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * Yโโป * ฯโป * ฮS - ฯโป * diffฮp * Yโโป * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ * Yโโบ * ฯโบ * ฮS + ฯโบ * diffฮp * Yโโบ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[1] * Yโโบ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[1] * Yโโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[1] * Yโโป * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[1] * Yโโบ * ฮS ))

        # v'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( ฯโบ * CโUโ * face.nฬ[2] * Yโโบ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโป * CโUแตฃ * face.nฬ[2] * Yโโป * ฮS ))
        
        A_vals[iแตฃ] -= ( ฯโป * CโUแตฃ * face.nฬ[2] * Yโโป * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโบ * CโUโ * face.nฬ[2] * Yโโบ * ฮS ))

        
        # T'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโTโ * Yโโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ * Yโโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * Yโโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ * Yโโบ * ฯโบ * ฮS ))

        
        # Yโ'
        iโ += 1; iแตฃ += 1

        A_vals[iโ] += ( Wโ * โฯโYโโ * Yโโบ * ฯโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * Yโโป * ฯโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS ))
        
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * Yโโป * ฯโป * ฮS + Wแตฃ * ฯโป * ฯโป * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * Yโโบ * ฯโบ * ฮS + Wโ * ฯโบ * ฯโบ * ฮS ))
        


        # ----------------------------
        flux = ( ฯโบ * Wโ * ฯโ + ฯโป * Wแตฃ * ฯแตฃ ) * ฮS
        B[ijStartโ + 1] -= flux
        B[ijStartแตฃ + 1] += flux

        flux = ( ฯโบ * Wโ * ฯโ * uโ + ฯโป * Wแตฃ * ฯแตฃ * uแตฃ ) * ฮS + ( pโ ) * face.nฬ[1] * ฮS
        B[ijStartโ + 2] -= flux
        B[ijStartแตฃ + 2] += flux

        flux = ( ฯโบ * Wโ * ฯโ * vโ + ฯโป * Wแตฃ * ฯแตฃ * vแตฃ ) * ฮS + ( pโ ) * face.nฬ[2] * ฮS
        B[ijStartโ + 3] -= flux
        B[ijStartแตฃ + 3] += flux

        flux = ( ฯโบ * Wโ * ฯโ * Hโโ + ฯโป * Wแตฃ * ฯแตฃ * Hโแตฃ ) * ฮS 
        B[ijStartโ + 4] -= flux
        B[ijStartแตฃ + 4] += flux

        flux = ( ฯโบ * Wโ * ฯโ * Yโโ + ฯโป * Wแตฃ * ฯแตฃ * Yโแตฃ ) * ฮS 
        B[ijStartโ + 5] -= flux
        B[ijStartแตฃ + 5] += flux




    end

    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    #bc_wall = []
    #append!( bc_wall, faces_boundary_top )
    #append!( bc_wall, faces_boundary_bottom )
    #append!( bc_wall, faces_boundary_left )
    #append!( bc_wall, faces_boundary_right )


    coupled_boundary!(
        ๐,cells,faces,
        faces_boundary_top, 
        ๐.top_p_BCtype, ๐.top_p_BCValue, 
        ๐.top_u_BCtype, ๐.top_u_BCValue, 
        ๐.top_v_BCtype, ๐.top_v_BCValue, 
        ๐.top_T_BCtype, ๐.top_T_BCValue, 
        ๐.top_Y_BCtype, ๐.top_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_bottom, 
        ๐.bottom_p_BCtype, ๐.bottom_p_BCValue, 
        ๐.bottom_u_BCtype, ๐.bottom_u_BCValue, 
        ๐.bottom_v_BCtype, ๐.bottom_v_BCValue, 
        ๐.bottom_T_BCtype, ๐.bottom_T_BCValue, 
        ๐.bottom_Y_BCtype, ๐.bottom_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_left, 
        ๐.left_p_BCtype, ๐.left_p_BCValue, 
        ๐.left_u_BCtype, ๐.left_u_BCValue, 
        ๐.left_v_BCtype, ๐.left_v_BCValue, 
        ๐.left_T_BCtype, ๐.left_T_BCValue, 
        ๐.left_Y_BCtype, ๐.left_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_right, 
        ๐.right_p_BCtype, ๐.right_p_BCValue, 
        ๐.right_u_BCtype, ๐.right_u_BCValue, 
        ๐.right_v_BCtype, ๐.right_v_BCValue, 
        ๐.right_T_BCtype, ๐.right_T_BCValue, 
        ๐.right_Y_BCtype, ๐.right_Y_BCValue,
        B_n, A_n, A_vals, B)


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ฮQ = solve(ps, A, B)

    #ฮQ = A\B

    #ml = ruge_stuben(A)
    #ฮQ = solve(ml, A)
    #P = aspreconditioner(ml)
    #ฮQ = bicgstabl(A, B, Pl = P)
    #ฮQ = gmres(A, B)

    #ฮQ = A\B

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

       cell.var[๐.p] += relax_p * ฮQ[ijStart + 1]
       cell.var[๐.u] += relax_U * ฮQ[ijStart + 2]
       cell.var[๐.v] += relax_U * ฮQ[ijStart + 3]
       cell.var[๐.T] += relax_T * ฮQ[ijStart + 4]
       cell.var[๐.Yโ] += relax_Y * ฮQ[ijStart + 5]

       cell.var[๐.p] = max(cell.var[๐.p],1.e-200)
       cell.var[๐.T] = max(cell.var[๐.T],1.e-200)
       
       norm_total += cell.var[๐.p]^2
       norm_total += cell.var[๐.u]^2
       norm_total += cell.var[๐.v]^2
       norm_total += cell.var[๐.T]^2
       norm_total += cell.var[๐.Yโ]^2

       diagon += 1
   end

   norm_total = sqrt(norm_total)

   #sleep(1000.0)


   return norm(ฮQ)/length(cells)/5
  

end



function coupled_boundary!(
    ๐::controls,
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
        
        ijStartโ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ฯโ = cells[face.owner].var[๐.ฯ]
        โฯโpโ = cells[face.owner].var[๐.โฯโp]
        โฯโTโ = cells[face.owner].var[๐.โฯโT]
        โHโโpโ = cells[face.owner].var[๐.โHโโp]
        โHโโTโ = cells[face.owner].var[๐.โHโโT]
        โHโโYโโ = cells[face.owner].var[๐.โHโโYโ]
        โฯโYโโ = cells[face.owner].var[๐.โฯโYโ]
        pโ = cells[face.owner].var[๐.p]
        Hโโ = cells[face.owner].var[๐.Hโ]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Tโ = cells[face.owner].var[๐.T]

        ฮS = face.ฮS

        Uโ = 0.0
        Uโ += cells[face.owner].var[๐.u]*face.nฬ[1]
        Uโ += cells[face.owner].var[๐.v]*face.nฬ[2]
        Uโ += cells[face.owner].var[๐.w]*face.nฬ[3]

        uโ = cells[face.owner].var[๐.u] - Uโ * face.nฬ[1]
        vโ = cells[face.owner].var[๐.v] - Uโ * face.nฬ[2]
        wโ = 0.0#cells[face.owner].var[๐.w] - Uโ * face.nฬ[3]

        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        id = []
        push!(id,i)
        push!(id,i+5)
        push!(id,i+10)
        push!(id,i+15)
        push!(id,i+20)

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pโ = cells[face.owner].var[๐.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pโ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pโ = p_BCValue(๐.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uโ = cells[face.owner].var[๐.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uโ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uโ = uโ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uโ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uโ = u_BCValue(๐.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vโ = cells[face.owner].var[๐.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vโ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vโ = vโ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vโ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vโ = v_BCValue(๐.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tโ = cells[face.owner].var[๐.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tโ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tโ = T_BCValue(๐.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yโโ = cells[face.owner].var[๐.Yโ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yโโ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yโโ = Y_BCValue(๐.time)
        end
        
        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        ฯโ, Hโโ, cโ = faceEOS!(๐,pโ,uโ,vโ,wโ,Tโ,Yโโ)
        #=
        # continuity
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * ฮS)

        
        # x-momentum
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * uโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * Uโ * ฮS + ฯโ * uโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * uโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * uโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * uโ * Uโ * ฮS)

        
        # y-momentum
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * vโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * vโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * Uโ * ฮS + ฯโ * vโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * vโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * vโ * Uโ * ฮS)


        # energy
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโpโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * uโ * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * vโ * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโTโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโYโโ * ฮS)


        # massfraction
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * Yโโ * ฮS)# + ฯโ * Hโโ * ๐.ฮt/ฯโ / ฮLR * ฮS
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * Yโโ * ฮS + ฯโ * Uโ * ฮS)
=#

        B[ijStartโ + 1] -= ( ฯโ * Uโ * ฮS )
        B[ijStartโ + 2] -= ( ฯโ * uโ * Uโ * ฮS + pโ * face.nฬ[1] * ฮS )
        B[ijStartโ + 3] -= ( ฯโ * vโ * Uโ * ฮS + pโ * face.nฬ[2] * ฮS )
        B[ijStartโ + 4] -= ( ฯโ * Hโโ * Uโ * ฮS )
        B[ijStartโ + 5] -= ( ฯโ * Yโโ * Uโ * ฮS )
        

    end
 


end












function coupled_Ap_boundary!(
    ๐::controls,
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

        ฯโ = cells[face.owner].var[๐.ฯ]
        pโ = cells[face.owner].var[๐.p]
        Hโโ = cells[face.owner].var[๐.Hโ]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Tโ = cells[face.owner].var[๐.T]

        ฮS = face.ฮS

        Uโ = 0.0
        Uโ += cells[face.owner].var[๐.u]*face.nฬ[1]
        Uโ += cells[face.owner].var[๐.v]*face.nฬ[2]
        Uโ += cells[face.owner].var[๐.w]*face.nฬ[3]

        uโ = cells[face.owner].var[๐.u] - Uโ * face.nฬ[1]
        vโ = cells[face.owner].var[๐.v] - Uโ * face.nฬ[2]
        wโ = 0.0#cells[face.owner].var[๐.w] - Uโ * face.nฬ[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pโ = cells[face.owner].var[๐.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pโ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pโ = p_BCValue(๐.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uโ = cells[face.owner].var[๐.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uโ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uโ = uโ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uโ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uโ = u_BCValue(๐.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vโ = cells[face.owner].var[๐.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vโ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vโ = vโ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vโ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vโ = v_BCValue(๐.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tโ = cells[face.owner].var[๐.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tโ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tโ = T_BCValue(๐.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yโโ = cells[face.owner].var[๐.Yโ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yโโ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yโโ = Y_BCValue(๐.time)
        end
        
        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        ฯโ, Hโโ, cโ = faceEOS!(๐,pโ,uโ,vโ,wโ,Tโ,Yโโ)
        

        flux = ฯโ * Uโ * ฮS
        Ap[face.owner] += flux / cells[face.owner].ฮฉ


    end
 


end



























function push_A_conv_diff!(
    A_rows::Array{Int64},
    A_cols::Array{Int64},
    A_vals::Array{Float64},
    AiL::Int64, iL::Int64, jL::Int64,
    AiR::Int64, iR::Int64, jR::Int64,
    convfluxโ::Float64, difffluxโ::Float64, 
    convfluxแตฃ::Float64, difffluxแตฃ::Float64
)
    A_vals[AiL] += ( convfluxโ + difffluxโ )
    push!(A_rows, iL)
    push!(A_cols, jL)
    push!(A_vals, convfluxแตฃ + difffluxแตฃ)
    
    A_vals[AiR] -= ( convfluxแตฃ + difffluxแตฃ )
    push!(A_rows, iR)
    push!(A_cols, jR)
    push!(A_vals, -( convfluxโ + difffluxโ ))
end


#=
function coupled!(
    ๐::controls,
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

        ฮฉ = cell.ฮฉ
        ฮt = ๐.ฮt
        u = cell.var[๐.u]
        v = cell.var[๐.v]
        ฯ = cell.var[๐.ฯ]
        Hโ = cell.var[๐.Hโ]
        p = cell.var[๐.p]
        Yโ = cell.var[๐.Yโ]
        โHโโp = cell.var[๐.โHโโp]
        โHโโT = cell.var[๐.โHโโT]
        โฯโp = cell.var[๐.โฯโp]
        โฯโT = cell.var[๐.โฯโT]
        โฯโYโ = cell.var[๐.โฯโYโ]
        โHโโYโ = cell.var[๐.โHโโYโ]
        ฯโฟ = cell.var[๐.ฯโฟ]
        uโฟ = cell.var[๐.uโฟ]
        vโฟ = cell.var[๐.vโฟ]
        Hโโฟ = cell.var[๐.Hโโฟ]
        pโฟ = cell.var[๐.pโฟ]
        Yโโฟ = cell.var[๐.Yโโฟ]
        ฯโฟโปยน = cell.var[๐.ฯโฟโปยน]
        uโฟโปยน = cell.var[๐.uโฟโปยน]
        vโฟโปยน = cell.var[๐.vโฟโปยน]
        Hโโฟโปยน = cell.var[๐.Hโโฟโปยน]
        pโฟโปยน = cell.var[๐.pโฟโปยน]
        Yโโฟโปยน = cell.var[๐.Yโโฟโปยน]

        #println(pโฟ,uโฟ,vโฟ,Hโโฟ)
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = โฯโp*ฮฉ/ฮt
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = โฯโT*ฮฉ/ฮt
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 5
        A_vals[i] = โฯโYโ*ฮฉ/ฮt

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = โฯโp*ฮฉ/ฮt * u - โฯโp*๐.gravity[1]*ฮฉ

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = ฯ*ฮฉ/ฮt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = โฯโT*ฮฉ/ฮt * u - โฯโT*๐.gravity[1]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 5
        A_vals[i] = โฯโYโ*ฮฉ/ฮt * u - โฯโYโ*๐.gravity[1]*ฮฉ


        # y-momentum
        #g = -9.8
        #g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = โฯโp*ฮฉ/ฮt * v - โฯโp*๐.gravity[2]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = ฯ*ฮฉ/ฮt

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = โฯโT*ฮฉ/ฮt * v - โฯโT*๐.gravity[2]*ฮฉ
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 5
        A_vals[i] = โฯโYโ*ฮฉ/ฮt * v - โฯโYโ*๐.gravity[2]*ฮฉ




        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = โฯโp*ฮฉ/ฮt * Hโ + โHโโp*ฮฉ/ฮt * ฯ - ฮฉ/ฮt

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = u*ฮฉ/ฮt * ฯ
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = v*ฮฉ/ฮt * ฯ
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = โฯโT*ฮฉ/ฮt * Hโ + โHโโT*ฮฉ/ฮt * ฯ
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 5
        A_vals[i] = โฯโYโ*ฮฉ/ฮt * Hโ + โHโโYโ*ฮฉ/ฮt * ฯ



        # mass fraction
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 1
        A_vals[i] = โฯโp*ฮฉ/ฮt * Yโ

        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 4
        A_vals[i] = โฯโT*ฮฉ/ฮt * Yโ 
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 5
        A_vals[i] = โฯโYโ*ฮฉ/ฮt * Yโ + ฮฉ/ฮt * ฯ
        
        

        # B
        if ๐.temporal_discretizationScheme == "1st"
            B[ijStart + 1] = -(ฯ - ฯโฟ)*ฮฉ/ฮt
            B[ijStart + 2] = -(ฯ*u - ฯโฟ*uโฟ)*ฮฉ/ฮt + ฯ*๐.gravity[1]*ฮฉ 
            B[ijStart + 3] = -(ฯ*v - ฯโฟ*vโฟ)*ฮฉ/ฮt + ฯ*๐.gravity[2]*ฮฉ 
            B[ijStart + 4] = -(ฯ*Hโ - ฯโฟ*Hโโฟ)*ฮฉ/ฮt + (p - pโฟ)*ฮฉ/ฮt
            B[ijStart + 5] = -(ฯ*Yโ - ฯโฟ*Yโโฟ)*ฮฉ/ฮt
        elseif ๐.temporal_discretizationScheme == "2nd"
            B[ijStart + 1] = -(1.5*ฯ - 2.0*ฯโฟ + 0.5*ฯโฟโปยน)*ฮฉ/ฮt
            B[ijStart + 2] = -(1.5*ฯ*u - 2.0*ฯโฟ*uโฟ + 0.5*ฯโฟโปยน*uโฟโปยน)*ฮฉ/ฮt + ฯ*๐.gravity[1]*ฮฉ 
            B[ijStart + 3] = -(1.5*ฯ*v - 2.0*ฯโฟ*vโฟ + 0.5*ฯโฟโปยน*vโฟโปยน)*ฮฉ/ฮt + ฯ*๐.gravity[2]*ฮฉ 
            B[ijStart + 4] = -(1.5*ฯ*Hโ - 2.0*ฯโฟ*Hโโฟ + 0.5*ฯโฟโปยน*Hโโฟโปยน)*ฮฉ/ฮt + (1.5*p - 2.0*pโฟ + 0.5*pโฟโปยน)*ฮฉ/ฮt
            B[ijStart + 5] = -(1.5*ฯ*Yโ - 2.0*ฯโฟ*Yโโฟ + 0.5*ฯโฟโปยน*Yโโฟโปยน)*ฮฉ/ฮt
        end



        diagon += 1


    end

    
    โฮpโx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pโ = 0.5 * (cells[face.owner].var[๐.p] + cells[face.neighbour].var[๐.p])
        โฮpโx0[face.owner, 1] += pโ * face.nฬ[1] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 2] += pโ * face.nฬ[2] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 3] += pโ * face.nฬ[3] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.neighbour, 1] -= pโ * face.nฬ[1] * face.ฮS / cells[face.neighbour].ฮฉ
        โฮpโx0[face.neighbour, 2] -= pโ * face.nฬ[2] * face.ฮS / cells[face.neighbour].ฮฉ
        โฮpโx0[face.neighbour, 3] -= pโ * face.nฬ[3] * face.ฮS / cells[face.neighbour].ฮฉ
    end

    for face in faces_boundary
        pโ = cells[face.owner].var[๐.p]
        โฮpโx0[face.owner, 1] += pโ * face.nฬ[1] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 2] += pโ * face.nฬ[2] * face.ฮS / cells[face.owner].ฮฉ
        โฮpโx0[face.owner, 3] += pโ * face.nฬ[3] * face.ฮS / cells[face.owner].ฮฉ
    end






    
    Ap = zeros(Float64, length(cells))
    for face in faces_internal


        ฯโ = cells[face.owner].var[๐.ฯ]
        ฯแตฃ = cells[face.neighbour].var[๐.ฯ]
        pโ = cells[face.owner].var[๐.p]
        pแตฃ = cells[face.neighbour].var[๐.p]
        #uโ = cells[face.owner].var[๐.u]
        #uแตฃ = cells[face.neighbour].var[๐.u]
        #vโ = cells[face.owner].var[๐.v]
        #vแตฃ = cells[face.neighbour].var[๐.v]
        uโ = cells[face.owner].var[๐.u]
        uแตฃ = cells[face.neighbour].var[๐.u]
        vโ = cells[face.owner].var[๐.v]
        vแตฃ = cells[face.neighbour].var[๐.v]

        Uโโ = uโ * face.nฬ[1] + vโ * face.nฬ[2]
        Uโแตฃ = uแตฃ * face.nฬ[1] + vแตฃ * face.nฬ[2]
        Uโ = 0.5 * (Uโโ + Uโแตฃ)

        centerโ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerแตฃ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ฮLR = norm(centerแตฃ - centerโ)

        ฯหข = 1.0 / (0.5/ฯโ + 0.5/ฯแตฃ)
        dฬ = ๐.ฮt / ฯหข
        
        Wโ = 0.0
        Wแตฃ = 0.0
        #if ๐.spatial_discretizationScheme == "upwind"
            Wโ = 0.5 * (1.0 + sign(Uโ))
            Wแตฃ = 1.0 - Wโ
        #elseif ๐.spatial_discretizationScheme == "central"
        #    Wโ = 0.5
        #    Wแตฃ = 1.0 - Wโ
        #end

        ฯโ = Wโ * ฯโ + Wแตฃ * ฯแตฃ
        
        # Rhie-Chow
        #=
        Uโ += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 1] * face.nฬ[1]
        Uโ += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 2] * face.nฬ[2]
        Uโ += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 3] * face.nฬ[3]
        Uโ += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 1] * face.nฬ[1]
        Uโ += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 2] * face.nฬ[2]
        Uโ += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 3] * face.nฬ[3]
        Uโ -= dฬ * (pแตฃ-pโ) / ฮLR
        =#

        flux = ฯโ * Uโ * face.ฮS
        Ap[face.owner] += Wโ * flux / cells[face.owner].ฮฉ
        Ap[face.neighbour] -= Wแตฃ * flux / cells[face.neighbour].ฮฉ
    end

    coupled_Ap_boundary!(
    ๐,cells,faces,
    faces_boundary_top, 
    ๐.top_p_BCtype, ๐.top_p_BCValue, 
    ๐.top_u_BCtype, ๐.top_u_BCValue, 
    ๐.top_v_BCtype, ๐.top_v_BCValue, 
    ๐.top_T_BCtype, ๐.top_T_BCValue, 
    ๐.top_Y_BCtype, ๐.top_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(๐,cells,faces,
    faces_boundary_bottom, 
    ๐.bottom_p_BCtype, ๐.bottom_p_BCValue, 
    ๐.bottom_u_BCtype, ๐.bottom_u_BCValue, 
    ๐.bottom_v_BCtype, ๐.bottom_v_BCValue, 
    ๐.bottom_T_BCtype, ๐.bottom_T_BCValue, 
    ๐.bottom_Y_BCtype, ๐.bottom_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(๐,cells,faces,
    faces_boundary_left, 
    ๐.left_p_BCtype, ๐.left_p_BCValue, 
    ๐.left_u_BCtype, ๐.left_u_BCValue, 
    ๐.left_v_BCtype, ๐.left_v_BCValue, 
    ๐.left_T_BCtype, ๐.left_T_BCValue, 
    ๐.left_Y_BCtype, ๐.left_Y_BCValue,
    Ap)

    coupled_Ap_boundary!(๐,cells,faces,
    faces_boundary_right, 
    ๐.right_p_BCtype, ๐.right_p_BCValue, 
    ๐.right_u_BCtype, ๐.right_u_BCValue, 
    ๐.right_v_BCtype, ๐.right_v_BCValue, 
    ๐.right_T_BCtype, ๐.right_T_BCValue, 
    ๐.right_Y_BCtype, ๐.right_Y_BCValue,
    Ap)
    






    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartโ = B_n*(face.owner-1)
        ijStartแตฃ = B_n*(face.neighbour-1)

       # ฯโ = cells[face.owner].var[๐.ฯ]
       # ฯแตฃ = cells[face.neighbour].var[๐.ฯ]

        ฯโ = face.varโ[๐.ฯ]
        ฯแตฃ = face.varแตฃ[๐.ฯ]

        pO = cells[face.owner].var[๐.p]
        pN = cells[face.neighbour].var[๐.p]
        pโ = face.varโ[๐.p]
        pแตฃ = face.varแตฃ[๐.p]

       # uโ = cells[face.owner].var[๐.u]
       # uแตฃ = cells[face.neighbour].var[๐.u]
       # vโ = cells[face.owner].var[๐.v]
       # vแตฃ = cells[face.neighbour].var[๐.v]

        uโ = face.varโ[๐.u]
        uแตฃ = face.varแตฃ[๐.u]
        vโ = face.varโ[๐.v]
        vแตฃ = face.varแตฃ[๐.v]

        wโ = 0.0#cells[face.owner].var[๐.w]
        wแตฃ = 0.0#cells[face.neighbour].var[๐.w]
       #= Hโโ = cells[face.owner].var[๐.Hโ]
        Hโแตฃ = cells[face.neighbour].var[๐.Hโ]
        #ฮผโ = cells[face.owner].var[๐.ฮผ]
        #ฮผแตฃ = cells[face.neighbour].var[๐.ฮผ]
        โฯโpโ = cells[face.owner].var[๐.โฯโp]
        โฯโpแตฃ = cells[face.neighbour].var[๐.โฯโp]
        โฯโTโ = cells[face.owner].var[๐.โฯโT]
        โฯโTแตฃ = cells[face.neighbour].var[๐.โฯโT]
        โHโโpโ = cells[face.owner].var[๐.โHโโp]
        โHโโpแตฃ = cells[face.neighbour].var[๐.โHโโp]
        โHโโTโ = cells[face.owner].var[๐.โHโโT]
        โHโโTแตฃ = cells[face.neighbour].var[๐.โHโโT]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Yโแตฃ = cells[face.neighbour].var[๐.Yโ]
        โฯโYโโ = cells[face.owner].var[๐.โฯโYโ]
        โฯโYโแตฃ = cells[face.neighbour].var[๐.โฯโYโ]
        โHโโYโโ = cells[face.owner].var[๐.โHโโYโ]
        โHโโYโแตฃ = cells[face.neighbour].var[๐.โHโโYโ]
        =#
        
        Hโโ = face.varโ[๐.Hโ]
        Hโแตฃ = face.varแตฃ[๐.Hโ]
        โฯโpโ = face.varโ[๐.โฯโp]
        โฯโpแตฃ = face.varแตฃ[๐.โฯโp]
        โฯโTโ = face.varโ[๐.โฯโT]
        โฯโTแตฃ = face.varแตฃ[๐.โฯโT]
        โHโโpโ = face.varโ[๐.โHโโp]
        โHโโpแตฃ = face.varแตฃ[๐.โHโโp]
        โHโโTโ = face.varโ[๐.โHโโT]
        โHโโTแตฃ = face.varแตฃ[๐.โHโโT]
        Yโโ = face.varโ[๐.Yโ]
        Yโแตฃ = face.varแตฃ[๐.Yโ]
        โฯโYโโ = face.varโ[๐.โฯโYโ]
        โฯโYโแตฃ = face.varแตฃ[๐.โฯโYโ]
        โHโโYโโ = face.varโ[๐.โHโโYโ]
        โHโโYโแตฃ = face.varแตฃ[๐.โHโโYโ]
        cโ = face.varโ[๐.c]
        cแตฃ = face.varแตฃ[๐.c]

        Uโโ = uโ * face.nฬ[1] + vโ * face.nฬ[2]
        Uโแตฃ = uแตฃ * face.nฬ[1] + vแตฃ * face.nฬ[2]
        Uโ = 0.5 * (Uโโ + Uโแตฃ)
        ฮS = face.ฮS

        centerโ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerแตฃ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ฮLR = norm(centerแตฃ - centerโ)

        ฯหข = 1.0 / (0.5/ฯโ + 0.5/ฯแตฃ)
        #d = 0.5 * (1.0 / (Ap[face.owner]) + 1.0 / (Ap[face.neighbour]) )
        dฬ = ๐.ฮt / ฯหข
        #dฬ = d / (2.0 + ฯหข / ๐.ฮt * d)
        #if d>1.e9
        #    dฬ = ๐.ฮt / ฯหข
        #end
        
        # Rhie-Chow
        Uโ_RC = 0.0
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 1] * face.nฬ[1]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 2] * face.nฬ[2]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯโ * โฮpโx0[face.owner, 3] * face.nฬ[3]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 1] * face.nฬ[1]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 2] * face.nฬ[2]
        Uโ_RC += dฬ * ฯหข * 0.5 / ฯแตฃ * โฮpโx0[face.neighbour, 3] * face.nฬ[3]
        Uโ_RC -= dฬ * (pN-pO) / ฮLR

        RCdiffฮp = dฬ / ฮLR

        # before step
        ฯโโฟ = cells[face.owner].var[๐.ฯโฟ]
        ฯแตฃโฟ = cells[face.neighbour].var[๐.ฯโฟ]
        ฯหขโฟ = 1.0 / (0.5/ฯโโฟ + 0.5/ฯแตฃโฟ)
        Uโโโฟ = cells[face.owner].var[๐.uโฟ] * face.nฬ[1] + cells[face.owner].var[๐.vโฟ] * face.nฬ[2]
        Uโแตฃโฟ = cells[face.neighbour].var[๐.uโฟ] * face.nฬ[1] + cells[face.neighbour].var[๐.vโฟ] * face.nฬ[2]
        #Uโ_RC += dฬ * ฯหขโฟ / ๐.ฮt * ( face.Uโโฟ - 0.5 * (Uโโโฟ + Uโแตฃโฟ) )
        #Uโ += ( face.Uโโฟ - 0.5 * (Uโโโฟ + Uโแตฃโฟ) )

        # YYL riemann
        cฬ = 0.5*(cโ + cแตฃ)
        Mโ = Uโโ/cฬ
        Mแตฃ = Uโแตฃ/cฬ
        # calculate M+ and P+ for left state
        Mโโบ = M_func(Mโ,1.0,0.125)
        pโบ = pre_func(Mโ,1.0,0.1875)
        # calculate M- and P- for left state
        Mแตฃโป = M_func(Mแตฃ,-1.0,0.125)
        pโป = pre_func(Mแตฃ,-1.0,0.1875)
        KLR = sqrt(0.5*(uโ^2+vโ^2+wโ^2+uแตฃ^2+vแตฃ^2+wแตฃ^2))
        Mdash = min(1.0,KLR/cฬ)

        Cdiffฮp = 0.5*(1.0-Mdash)^2/cฬ /(0.5*(ฯโ+ฯแตฃ))
        WUโ = Mโโบ/Mโ
        if abs(Mโ) <= 0.00001
            WUโ = 0.5
        end
        WUแตฃ = Mแตฃโป/Mแตฃ
        if abs(Mแตฃ) <= 0.00001
            WUแตฃ = 0.5
        end

	    #Uโ = WUโ*Uโโ + WUแตฃ*Uโแตฃ + Uโ_RC - Cdiffฮp * (pแตฃ-pโ)
	    Uโ = WUโ*Uโโ + WUแตฃ*Uโแตฃ - Cdiffฮp * (pแตฃ-pโ)

        Wpโ = pโบ
        Wpแตฃ = pโป

        #diffฮp = RCdiffฮp + Cdiffฮp
        diffฮp = Cdiffฮp

        #--------------------
        # SAVE
        face.Uโ = Uโ
        #--------------------




        
        Wโ = 0.5 * (1.0 + sign(Uโ))
        Wแตฃ = 1.0 - Wโ

        ฯโ = Wโ * ฯโ + Wแตฃ * ฯแตฃ
        uโ = Wโ * uโ + Wแตฃ * uแตฃ
        vโ = Wโ * vโ + Wแตฃ * vแตฃ
        wโ = 0.0#Wโ * wโ + Wแตฃ * wแตฃ
        Hโโ = Wโ * Hโโ + Wแตฃ * Hโแตฃ
        Yโโ = Wโ * Yโโ + Wแตฃ * Yโแตฃ

        pโ = Wpโ*pโ + Wpแตฃ*pแตฃ

        
        iโ = A_n*(face.owner-1)
        iแตฃ = A_n*(face.neighbour-1)


        #--- ACID ----
        ฯโ_ACID, ฯแตฃ_ACID, โฯโpโ_ACID, โฯโpแตฃ_ACID, โHโโpโ_ACID, โHโโpแตฃ_ACID,
        Hโโ_ACID, Hโแตฃ_ACID, โฯโTโ_ACID, โฯโTแตฃ_ACID, โHโโTโ_ACID, โHโโTแตฃ_ACID,
        โฯโYโโ_ACID, โฯโYโแตฃ_ACID =
        EOS_ACID(
            ๐,
            pโ,pแตฃ,
            cells[face.owner].var[๐.u],cells[face.neighbour].var[๐.u],
            cells[face.owner].var[๐.v],cells[face.neighbour].var[๐.v],
            cells[face.owner].var[๐.w],cells[face.neighbour].var[๐.w],
            cells[face.owner].var[๐.T],cells[face.neighbour].var[๐.T],
            cells[face.owner].var[๐.ฮฑโ],cells[face.neighbour].var[๐.ฮฑโ]
            
        )

        โฯโpโ_ACID = โฯโpโ
        โฯโpแตฃ_ACID = โฯโpแตฃ
        โHโโpโ_ACID = โHโโpโ
        โHโโpแตฃ_ACID = โHโโpแตฃ
        โฯโTโ_ACID = โฯโTโ
        โฯโTแตฃ_ACID = โฯโTแตฃ
        โHโโTโ_ACID = โHโโTโ
        โHโโTแตฃ_ACID = โHโโTแตฃ

        ฯโ_ACID = ฯโ
        ฯแตฃ_ACID = ฯแตฃ
        Hโโ_ACID = Hโโ
        Hโแตฃ_ACID = Hโแตฃ
        #=
        =#

        #=
        ฯโโ_ACID = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        ฯโแตฃ_ACID = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        
        Hโโโ_ACID = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        Hโโแตฃ_ACID = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        =#



        #------------------------
        # continuity
        # p'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโpโ * Uโ * ฮS + ฯโ * diffฮp * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ_ACID * Uโ * ฮS - ฯโ * diffฮp * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * Uโ * ฮS - ฯโ * diffฮp * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ_ACID * Uโ * ฮS + ฯโ * diffฮp * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1
        
        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[1] * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[1] * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[1] * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[1] * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[2] * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[2] * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[2] * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[2] * ฮS ))
        
        # T'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโTโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ_ACID * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ_ACID * Uโ * ฮS ))
        
        # Yโ'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโYโโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 1); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 1); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * Uโ * ฮS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโpโ * uโ * Uโ * ฮS + Wpโ * face.nฬ[1] * ฮS + ฯโ * diffฮp * uโ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ_ACID * uโ * Uโ * ฮS + Wpแตฃ * face.nฬ[1] * ฮS - ฯโ * diffฮp * uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * uโ * Uโ * ฮS + Wpแตฃ * face.nฬ[1] * ฮS - ฯโ * diffฮp * uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ_ACID * uโ * Uโ * ฮS + Wpโ * face.nฬ[1] * ฮS + ฯโ * diffฮp * uโ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[1] * uโ * ฮS + Wโ * ฯโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[1] * uโ * ฮS + Wแตฃ * ฯแตฃ_ACID * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[1] * uโ * ฮS + Wแตฃ * ฯแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[1] * uโ * ฮS + Wโ * ฯโ_ACID * Uโ * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[2] * uโ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[2] * uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[2] * uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[2] * uโ * ฮS ))

        # T'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโTโ * uโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ_ACID * uโ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * uโ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ_ACID * uโ * Uโ * ฮS ))

        # Yโ'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโYโโ * uโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 2); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * uโ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * uโ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 2); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * uโ * Uโ * ฮS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] +=  ( Wโ * โฯโpโ* vโ * Uโ * ฮS + Wpโ * face.nฬ[2] * ฮS + ฯโ * diffฮp * vโ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ_ACID * vโ * Uโ * ฮS + Wpแตฃ * face.nฬ[2] * ฮS - ฯโ * diffฮp * vโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ* vโ * Uโ * ฮS + Wpแตฃ * face.nฬ[2] * ฮS - ฯโ * diffฮp * vโ * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ_ACID * vโ * Uโ * ฮS + Wpโ * face.nฬ[2] * ฮS + ฯโ * diffฮp * vโ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[1] * vโ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[1] * vโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[1] * vโ * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[1] * vโ * ฮS ))
        
        # v'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[2] * vโ * ฮS + Wโ * ฯโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[2] * vโ * ฮS + Wแตฃ * ฯแตฃ_ACID * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[2] * vโ * ฮS + Wแตฃ * ฯแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[2] * vโ * ฮS + Wโ * ฯโ_ACID * Uโ * ฮS ))

        # T'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโTโ * vโ * Uโ *ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ_ACID * vโ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * vโ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -(  Wโ * โฯโTโ_ACID * vโ * Uโ *ฮS ))

        # Yโ'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโYโโ * vโ * Uโ *ฮS )
        push!(A_rows, ijStartโ + 3); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * vโ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * vโ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 3); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -(  Wโ * โฯโYโโ * vโ * Uโ *ฮS ))


        

        #------------------------
        # energy
        # p'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโpโ * Hโโ * Uโ * ฮS + ฯโ * Wโ * โHโโpโ * Uโ * ฮS + ฯโ * diffฮp * Hโโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ_ACID * Hโโ * Uโ * ฮS + ฯโ * Wแตฃ * โHโโpแตฃ_ACID * Uโ * ฮS - ฯโ * diffฮp * Hโโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * Hโโ * Uโ * ฮS + ฯโ * Wแตฃ * โHโโpแตฃ * Uโ * ฮS - ฯโ * diffฮp * Hโโ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ_ACID * Hโโ * Uโ * ฮS + ฯโ * Wโ * โHโโpโ_ACID * Uโ * ฮS + ฯโ * diffฮp * Hโโ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * Wโ * uโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * Wแตฃ * uแตฃ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * Wแตฃ * uแตฃ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * Wโ * uโ * ฮS ))

        # v'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * Wโ * vโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * Wแตฃ * vแตฃ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * Wแตฃ * vแตฃ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * Wโ * vโ * ฮS ))

        
        # T'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโTโ * Hโโ * Uโ * ฮS + ฯโ * Wโ * โHโโTโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ_ACID * Hโแตฃ_ACID * Uโ * ฮS + ฯโ * Wแตฃ * โHโโTแตฃ_ACID * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * Hโแตฃ * Uโ * ฮS + ฯโ * Wแตฃ * โHโโTแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ_ACID * Hโโ_ACID * Uโ * ฮS + ฯโ * Wโ * โHโโTโ_ACID * Uโ * ฮS ))

        
        # Yโ'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโYโโ * Hโโ * Uโ * ฮS + ฯโ * Wโ * โHโโYโโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 4); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * Hโแตฃ_ACID * Uโ * ฮS + ฯโ * Wแตฃ * โHโโYโแตฃ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * Hโแตฃ * Uโ * ฮS + ฯโ * Wแตฃ * โHโโYโแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 4); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * Hโโ_ACID * Uโ * ฮS + ฯโ * Wโ * โHโโYโโ * Uโ * ฮS ))
        

        

        #------------------------
        # mass fraction
        # p'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโpโ * Yโโ * Uโ * ฮS + ฯโ * diffฮp * Yโโ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 1)
        push!(A_vals, ( Wแตฃ * โฯโpแตฃ_ACID * Yโโ * Uโ * ฮS - ฯโ * diffฮp * Yโโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโpแตฃ * Yโโ * Uโ * ฮS - ฯโ * diffฮp * Yโโ * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 1)
        push!(A_vals, -( Wโ * โฯโpโ_ACID * Yโโ * Uโ * ฮS + ฯโ * diffฮp * Yโโ * ฮS ))
        
        # u'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[1] * Yโโ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 2)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[1] * Yโโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[1] * Yโโ * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 2)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[1] * Yโโ * ฮS ))

        # v'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        A_vals[iโ] += ( ฯโ * WUโ * face.nฬ[2] * Yโโ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 3)
        push!(A_vals, ( ฯโ * WUแตฃ * face.nฬ[2] * Yโโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        A_vals[iแตฃ] -= ( ฯโ * WUแตฃ * face.nฬ[2] * Yโโ * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 3)
        push!(A_vals, -( ฯโ * WUโ * face.nฬ[2] * Yโโ * ฮS ))

        
        # T'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโTโ * Yโโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 4)
        push!(A_vals, ( Wแตฃ * โฯโTแตฃ_ACID * Yโโ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโTแตฃ * Yโโ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 4)
        push!(A_vals, -( Wโ * โฯโTโ_ACID * Yโโ * Uโ * ฮS ))

        
        # Yโ'
        iโ += 1; iแตฃ += 1

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        A_vals[iโ] += ( Wโ * โฯโYโโ * Yโโ * Uโ * ฮS + ฯโ * Wโ * Uโ * ฮS )
        push!(A_rows, ijStartโ + 5); push!(A_cols, ijStartแตฃ + 5)
        push!(A_vals, ( Wแตฃ * โฯโYโแตฃ * Yโโ * Uโ * ฮS + ฯโ * Wแตฃ * Uโ * ฮS ))
        
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        A_vals[iแตฃ] -= ( Wแตฃ * โฯโYโแตฃ * Yโโ * Uโ * ฮS + ฯโ * Wแตฃ * Uโ * ฮS )
        push!(A_rows, ijStartแตฃ + 5); push!(A_cols, ijStartโ + 5)
        push!(A_vals, -( Wโ * โฯโYโโ * Yโโ * Uโ * ฮS + ฯโ * Wโ * Uโ * ฮS ))
        

        # ----------------------------

        # B

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        B[ijStartโ + 1] -= ( ฯโ * Uโ * ฮS )
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        B[ijStartแตฃ + 1] += ( ฯโ * Uโ * ฮS )
        # B

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        B[ijStartโ + 2] -= ( ฯโ * uโ * Uโ * ฮS + pโ * face.nฬ[1] * ฮS )
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        B[ijStartแตฃ + 2] += ( ฯโ * uโ * Uโ * ฮS + pโ * face.nฬ[1] * ฮS )
        # B

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        B[ijStartโ + 3] -= ( ฯโ * vโ * Uโ * ฮS + pโ * face.nฬ[2] * ฮS )
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        B[ijStartแตฃ + 3] += ( ฯโ * vโ * Uโ * ฮS + pโ * face.nฬ[2] * ฮS )
        # B

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        Hโโ = (Wโ * Hโโ + Wแตฃ * Hโแตฃ_ACID)
        B[ijStartโ + 4] -= ( ฯโ * Hโโ * Uโ * ฮS )
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        Hโโ = (Wโ * Hโโ_ACID + Wแตฃ * Hโแตฃ)
        B[ijStartแตฃ + 4] += ( ฯโ * Hโโ * Uโ * ฮS )
        # B

        ฯโ = (Wโ*ฯโ + Wแตฃ*ฯแตฃ_ACID)
        B[ijStartโ + 5] -= ( ฯโ * Yโโ * Uโ * ฮS )
        ฯโ = (Wโ*ฯโ_ACID + Wแตฃ*ฯแตฃ)
        B[ijStartแตฃ + 5] += ( ฯโ * Yโโ * Uโ * ฮS )




    end


    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    #bc_wall = []
    #append!( bc_wall, faces_boundary_top )
    #append!( bc_wall, faces_boundary_bottom )
    #append!( bc_wall, faces_boundary_left )
    #append!( bc_wall, faces_boundary_right )


    coupled_boundary!(
        ๐,cells,faces,
        faces_boundary_top, 
        ๐.top_p_BCtype, ๐.top_p_BCValue, 
        ๐.top_u_BCtype, ๐.top_u_BCValue, 
        ๐.top_v_BCtype, ๐.top_v_BCValue, 
        ๐.top_T_BCtype, ๐.top_T_BCValue, 
        ๐.top_Y_BCtype, ๐.top_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_bottom, 
        ๐.bottom_p_BCtype, ๐.bottom_p_BCValue, 
        ๐.bottom_u_BCtype, ๐.bottom_u_BCValue, 
        ๐.bottom_v_BCtype, ๐.bottom_v_BCValue, 
        ๐.bottom_T_BCtype, ๐.bottom_T_BCValue, 
        ๐.bottom_Y_BCtype, ๐.bottom_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_left, 
        ๐.left_p_BCtype, ๐.left_p_BCValue, 
        ๐.left_u_BCtype, ๐.left_u_BCValue, 
        ๐.left_v_BCtype, ๐.left_v_BCValue, 
        ๐.left_T_BCtype, ๐.left_T_BCValue, 
        ๐.left_Y_BCtype, ๐.left_Y_BCValue,
        B_n, A_n, A_vals, B)

    coupled_boundary!(๐,cells,faces,
        faces_boundary_right, 
        ๐.right_p_BCtype, ๐.right_p_BCValue, 
        ๐.right_u_BCtype, ๐.right_u_BCValue, 
        ๐.right_v_BCtype, ๐.right_v_BCValue, 
        ๐.right_T_BCtype, ๐.right_T_BCValue, 
        ๐.right_Y_BCtype, ๐.right_Y_BCValue,
        B_n, A_n, A_vals, B)


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ฮQ = solve(ps, A, B)

    #ฮQ = A\B

    #ml = ruge_stuben(A)
    #ฮQ = solve(ml, A)
    #P = aspreconditioner(ml)
    #ฮQ = bicgstabl(A, B, Pl = P)
    #ฮQ = gmres(A, B)

    #ฮQ = A\B

   # error()
   # exit()

    relax_p = 0.9
    relax_U = 0.9
    relax_T = 0.9
    relax_Y = 0.9


    diagon = 1
    maximum_p = -1.e12
    maximum_U = -1.e12
    maximum_T = -1.e12
    maximum_Y = -1.e12
    norm_p = 0.0
    norm_U = 0.0
    norm_T = 0.0
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[๐.p] += relax_p * ฮQ[ijStart + 1]
        cell.var[๐.u] += relax_U * ฮQ[ijStart + 2]
        cell.var[๐.v] += relax_U * ฮQ[ijStart + 3]
        cell.var[๐.T] += relax_T * ฮQ[ijStart + 4]
        cell.var[๐.Yโ] += relax_Y * ฮQ[ijStart + 5]

        cell.var[๐.p] = max(cell.var[๐.p],1.e-200)
        cell.var[๐.T] = max(cell.var[๐.T],1.e-200)
        
        #println(cell.var[๐.p])
        
        norm_p += ฮQ[ijStart + 1]^2
        norm_U += ฮQ[ijStart + 2]^2
        norm_U += ฮQ[ijStart + 3]^2
        norm_T += ฮQ[ijStart + 4]^2
        maximum_p = max(maximum_p,abs(cell.var[๐.p]))
        maximum_U = max(maximum_U,abs(cell.var[๐.u]))
        maximum_U = max(maximum_U,abs(cell.var[๐.v]))
        maximum_T = max(maximum_T,abs(cell.var[๐.T]))
        maximum_Y = max(maximum_Y,abs(cell.var[๐.Yโ]))

        diagon += 1
    end

    #sleep(1000.0)


    return norm(ฮQ),log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-200)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-200)),
    log10(sqrt(norm_T)/length(cells)/(maximum_T+1.e-200))
   


end



function coupled_boundary!(
    ๐::controls,
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
        
        ijStartโ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ฯโ = cells[face.owner].var[๐.ฯ]
        โฯโpโ = cells[face.owner].var[๐.โฯโp]
        โฯโTโ = cells[face.owner].var[๐.โฯโT]
        โHโโpโ = cells[face.owner].var[๐.โHโโp]
        โHโโTโ = cells[face.owner].var[๐.โHโโT]
        โHโโYโโ = cells[face.owner].var[๐.โHโโYโ]
        โฯโYโโ = cells[face.owner].var[๐.โฯโYโ]
        pโ = cells[face.owner].var[๐.p]
        Hโโ = cells[face.owner].var[๐.Hโ]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Tโ = cells[face.owner].var[๐.T]

        ฮS = face.ฮS

        Uโ = 0.0
        Uโ += cells[face.owner].var[๐.u]*face.nฬ[1]
        Uโ += cells[face.owner].var[๐.v]*face.nฬ[2]
        Uโ += cells[face.owner].var[๐.w]*face.nฬ[3]

        uโ = cells[face.owner].var[๐.u] - Uโ * face.nฬ[1]
        vโ = cells[face.owner].var[๐.v] - Uโ * face.nฬ[2]
        wโ = 0.0#cells[face.owner].var[๐.w] - Uโ * face.nฬ[3]

        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        id = []
        push!(id,i)
        push!(id,i+5)
        push!(id,i+10)
        push!(id,i+15)
        push!(id,i+20)

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pโ = cells[face.owner].var[๐.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pโ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pโ = p_BCValue(๐.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uโ = cells[face.owner].var[๐.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uโ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uโ = uโ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uโ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uโ = u_BCValue(๐.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vโ = cells[face.owner].var[๐.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vโ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vโ = vโ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vโ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vโ = v_BCValue(๐.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tโ = cells[face.owner].var[๐.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tโ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tโ = T_BCValue(๐.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yโโ = cells[face.owner].var[๐.Yโ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yโโ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yโโ = Y_BCValue(๐.time)
        end
        
        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        ฯโ, Hโโ, cโ = faceEOS!(๐,pโ,uโ,vโ,wโ,Tโ,Yโโ)
        
        # continuity
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * ฮS)

        
        # x-momentum
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * uโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * Uโ * ฮS + ฯโ * uโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * uโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * uโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * uโ * Uโ * ฮS)

        
        # y-momentum
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * vโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * vโ * face.nฬ[1] * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * Uโ * ฮS + ฯโ * vโ * face.nฬ[2] * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * vโ * Uโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * vโ * Uโ * ฮS)


        # energy
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโpโ * ฮS)
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * Hโโ * ฮS + ฯโ * Uโ * uโ * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * Hโโ * ฮS + ฯโ * Uโ * vโ * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโTโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * Hโโ * ฮS + ฯโ * Uโ * โHโโYโโ * ฮS)


        # massfraction
        i += 1
        A_vals[i] += coeff_p * (โฯโpโ * Uโ * Yโโ * ฮS)# + ฯโ * Hโโ * ๐.ฮt/ฯโ / ฮLR * ฮS
        i += 1
        A_vals[i] += coeff_u * (ฯโ * face.nฬ[1] * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_v * (ฯโ * face.nฬ[2] * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_T * (โฯโTโ * Uโ * Yโโ * ฮS)
        i += 1
        A_vals[i] += coeff_Y * (โฯโYโโ * Uโ * Yโโ * ฮS + ฯโ * Uโ * ฮS)


        B[ijStartโ + 1] -= ( ฯโ * Uโ * ฮS )
        B[ijStartโ + 2] -= ( ฯโ * uโ * Uโ * ฮS + pโ * face.nฬ[1] * ฮS )
        B[ijStartโ + 3] -= ( ฯโ * vโ * Uโ * ฮS + pโ * face.nฬ[2] * ฮS )
        B[ijStartโ + 4] -= ( ฯโ * Hโโ * Uโ * ฮS )
        B[ijStartโ + 5] -= ( ฯโ * Yโโ * Uโ * ฮS )
        

    end
 


end












function coupled_Ap_boundary!(
    ๐::controls,
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

        ฯโ = cells[face.owner].var[๐.ฯ]
        pโ = cells[face.owner].var[๐.p]
        Hโโ = cells[face.owner].var[๐.Hโ]
        Yโโ = cells[face.owner].var[๐.Yโ]
        Tโ = cells[face.owner].var[๐.T]

        ฮS = face.ฮS

        Uโ = 0.0
        Uโ += cells[face.owner].var[๐.u]*face.nฬ[1]
        Uโ += cells[face.owner].var[๐.v]*face.nฬ[2]
        Uโ += cells[face.owner].var[๐.w]*face.nฬ[3]

        uโ = cells[face.owner].var[๐.u] - Uโ * face.nฬ[1]
        vโ = cells[face.owner].var[๐.v] - Uโ * face.nฬ[2]
        wโ = 0.0#cells[face.owner].var[๐.w] - Uโ * face.nฬ[3]

        coeff_p = 0.0
        if p_BCtype == "zeroGradient"
            coeff_p = 1.0
            pโ = cells[face.owner].var[๐.p]
        elseif p_BCtype == "fixedValue"
            coeff_p = 0.0
            pโ = p_BCValue
        elseif p_BCtype == "function"
            coeff_p = 0.0
            pโ = p_BCValue(๐.time)
        end
        
        coeff_u = 0.0
        if u_BCtype == "zeroGradient"
            coeff_u = 1.0
            uโ = cells[face.owner].var[๐.u]
        elseif u_BCtype == "fixedValue"
            coeff_u = 0.0
            uโ = u_BCValue
        elseif u_BCtype == "slip"
            coeff_u = 0.0
            uโ = uโ
        elseif u_BCtype == "wall"
            coeff_u = 0.0
            uโ = 0.0
        elseif u_BCtype == "function"
            coeff_u = 0.0
            uโ = u_BCValue(๐.time)
        end
        
        coeff_v = 0.0
        if v_BCtype == "zeroGradient"
            coeff_v = 1.0
            vโ = cells[face.owner].var[๐.v]
        elseif v_BCtype == "fixedValue"
            coeff_v = 0.0
            vโ = v_BCValue
        elseif v_BCtype == "slip"
            coeff_v = 0.0
            vโ = vโ
        elseif v_BCtype == "wall"
            coeff_v = 0.0
            vโ = 0.0
        elseif v_BCtype == "function"
            coeff_v = 0.0
            vโ = v_BCValue(๐.time)
        end
        
        coeff_T = 0.0
        if T_BCtype == "zeroGradient"
            coeff_T = 1.0
            Tโ = cells[face.owner].var[๐.T]
        elseif T_BCtype == "fixedValue"
            coeff_T = 0.0
            Tโ = T_BCValue
        elseif T_BCtype == "function"
            coeff_T = 0.0
            Tโ = T_BCValue(๐.time)
        end
        
        coeff_Y = 0.0
        if T_BCtype == "zeroGradient"
            coeff_Y = 1.0
            Yโโ = cells[face.owner].var[๐.Yโ]
        elseif Y_BCtype == "fixedValue"
            coeff_Y = 0.0
            Yโโ = Y_BCValue
        elseif Y_BCtype == "function"
            coeff_Y = 0.0
            Yโโ = Y_BCValue(๐.time)
        end
        
        Uโ = uโ * face.nฬ[1] + vโ * face.nฬ[2] + wโ * face.nฬ[3]

        ฯโ, Hโโ, cโ = faceEOS!(๐,pโ,uโ,vโ,wโ,Tโ,Yโโ)
        

        flux = ฯโ * Uโ * ฮS
        Ap[face.owner] += flux / cells[face.owner].ฮฉ


    end
 


end


=#


function M_func(M::Float64, op::Float64, ฮฑ::Float64)
    mu = 0.0
	if abs(M) > 1.0 
		mu = 0.5*(M + op*abs(M))
	else
		mu = op*0.25*(M + op)^2.0 + op*ฮฑ*(M*M-1.0)^2.0
    end
	
	return mu
end

function pre_func(M::Float64, op::Float64, ฮฑ::Float64)
    mu = 0.0
	if abs(M) > 1.0
		mu = 0.5*(1.0 + op*sign(M) )
	else
		mu = 0.25*(M + op)^2.0*(2.0-op*M) + op*ฮฑ*M*(M*M-1.0)^2.0
    end
	
	return mu;
end


