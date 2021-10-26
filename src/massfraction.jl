

function massfraction!(
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




    #A_rows::Vector{Int64} = []
    #A_cols::Vector{Int64} = []
    #A_vals::Vector{Float64} = []
    
    # contruct A matrix diagonal terms
    # contruct B vector
    A_rows = zeros(Int64, length(cells))
    A_cols = zeros(Int64, length(cells))
    A_vals = zeros(Float64, length(cells))
    B = zeros(Float64, length(cells))
    
    diagon = 1

    for cell in cells

        ρ = cell.var[👉.ρ]
        Y₁ = cell.var[👉.Y₁]
        ρⁿ = cell.var[👉.ρⁿ]
        Y₁ⁿ = cell.var[👉.Y₁ⁿ]
        ∂ρ∂Y₁ = cell.var[👉.∂ρ∂Y₁]
        ρⁿ⁻¹ = cell.var[👉.ρⁿ⁻¹]
        Y₁ⁿ⁻¹ = cell.var[👉.Y₁ⁿ⁻¹]
        Ω = cell.Ω
        Δt = 👉.Δt
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon

        A_vals[diagon] = ρ*Ω/Δt
        A_vals[diagon] += ∂ρ∂Y₁*Y₁*Ω/Δt

        if 👉.temporal_discretizationScheme == "1st"
            B[diagon] = -(ρ*Y₁ - ρⁿ*Y₁ⁿ)*Ω/Δt
        elseif 👉.temporal_discretizationScheme == "2nd"
            B[diagon] = -(1.5*ρ*Y₁ - 2.0*ρⁿ*Y₁ⁿ + 0.5*ρⁿ⁻¹*Y₁ⁿ⁻¹)*Ω/Δt
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







    ∂Y₁∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.Y₁] + cells[face.neighbour].var[👉.Y₁])
        ∂Y₁∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Y₁∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Y₁∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂Y₁∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂Y₁∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂Y₁∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.Y₁]
        ∂Y₁∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Y₁∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Y₁∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_Y == "upwind"
        calcUpwind!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "central"
        calcCentral!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "quick"
        calcQUICK!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "minmod"
        calcMinmod!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "vanleer"
        calcVanLeer!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "superbee"
        calcSUPERBEE!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "linear"
        calcLinearInterpolation!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    end




    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
#=
        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        ∂ρ∂Y₁ₗ = cells[face.owner].var[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂ρ∂Y₁]
=#
        pO = cells[face.owner].var[👉.p]
        pN = cells[face.neighbour].var[👉.p]

        ρₗ = face.varₗ[👉.ρ]
        ρᵣ = face.varᵣ[👉.ρ]
        pₗ = face.varₗ[👉.p]
        pᵣ = face.varᵣ[👉.p]
        uₗ = face.varₗ[👉.u]
        uᵣ = face.varᵣ[👉.u]
        vₗ = face.varₗ[👉.v]
        vᵣ = face.varᵣ[👉.v]
        ∂ρ∂Y₁ₗ = face.varₗ[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = face.varᵣ[👉.∂ρ∂Y₁]


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
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pN-pO) / ΔLR

        
        # before step
        ρₗⁿ = cells[face.owner].var[👉.ρⁿ]
        ρᵣⁿ = cells[face.neighbour].var[👉.ρⁿ]
        ρˢⁿ = 1.0 / (0.5/ρₗⁿ + 0.5/ρᵣⁿ)
        Uₙₗⁿ = cells[face.owner].var[👉.uⁿ] * face.n̂[1] + cells[face.owner].var[👉.vⁿ] * face.n̂[2]
        Uₙᵣⁿ = cells[face.neighbour].var[👉.uⁿ] * face.n̂[1] + cells[face.neighbour].var[👉.vⁿ] * face.n̂[2]
        Uₙ += d̂ * ρˢⁿ / 👉.Δt * ( face.Uₙⁿ - 0.5 * (Uₙₗⁿ + Uₙᵣⁿ) )
        #Uₙ += ( face.Uₙⁿ - 0.5 * (Uₙₗⁿ + Uₙᵣⁿ) )
        #--------------------
        # SAVE
        face.Uₙ = Uₙ
        #--------------------


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
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        
        Y₁ₗ = cells[face.owner].var[👉.Y₁]
        Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]
        Y₁ₙ = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ

        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( Wₗ * ρₙ * Uₙ * ΔS + Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS )
        push!(A_vals, ( Wᵣ * ρₙ * Uₙ * ΔS + Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS ))
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( Wᵣ * ρₙ * Uₙ * ΔS + Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS )
        push!(A_vals, -( Wₗ * ρₙ * Uₙ * ΔS + Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS ) )


        # convective terms
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        B[face.neighbour] += ρₙ * Y₁ₙ * Uₙ * ΔS


    end
    

    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    bc_wall = []
    append!( bc_wall, faces_boundary_top )
    append!( bc_wall, faces_boundary_bottom )
    append!( bc_wall, faces_boundary_left )
    append!( bc_wall, faces_boundary_right )

    bc_slipwall = []
    
    bc_subinlet = []
    #append!( bc_subinlet, faces_boundary_left )
    
    bc_suboutlet = []
    #append!( bc_supoutlet, faces_boundary_right )
    
    bc_supoutlet = []
    
    for face in bc_wall
        ΔS = face.ΔS

        uₙ = 0.0
        vₙ = 0.0
        wₙ = 0.0
        Uₙ = 0.0
        ρₙ = cells[face.owner].var[👉.ρ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        A_vals[face.owner] += ρₙ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        
    end
    

    for face in bc_slipwall
        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]

        ρₙ = cells[face.owner].var[👉.ρ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        A_vals[face.owner] += ρₙ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        
    end
    
    for face in bc_subinlet
        ΔS = face.ΔS

        uₙ = 1.0
        vₙ = 0.0

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        ρₙ = cells[face.owner].var[👉.ρ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        #A_vals[face.owner] += ρₙ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        
    end

    for face in bc_suboutlet
        ΔS = face.ΔS


        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        ρₙ = cells[face.owner].var[👉.ρ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        A_vals[face.owner] += ρₙ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        
    end
    

    for face in bc_supoutlet
        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        ρₙ = cells[face.owner].var[👉.ρ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        A_vals[face.owner] += ρₙ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= ρₙ * Y₁ₙ * Uₙ * ΔS
        
    end


 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ΔY₁ = solve(ps, A, B)
    





    relax = 0.9




    diagon = 1
    maximum_Y = -1.e12
    norm_Y = 0.0
    for cell in cells

        cell.var[👉.Y₁] += relax*ΔY₁[diagon]

        norm_Y += ΔY₁[diagon]^2
        maximum_Y = max(maximum_Y,abs(cell.var[👉.Y₁]))

        diagon += 1
    end


    return norm(ΔY₁)
   
   

end
