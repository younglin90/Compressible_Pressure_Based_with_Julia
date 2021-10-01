

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
        Ω = cell.Ω
        Δt = 👉.Δt
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon

        A_vals[diagon] = ρ*Ω/Δt
        A_vals[diagon] += ∂ρ∂Y₁*Y₁*Ω/Δt

        B[diagon] = -(ρ*Y₁ - ρⁿ*Y₁ⁿ)*Ω/Δt

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

    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal

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
        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        #invρΔt = (wₗ/ρₗ + wᵣ/ρᵣ) * 👉.Δt
        invρΔt = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        #Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        #Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        #=
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        =#
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ
        
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

#=
        #--- ACID ----
        ρₗ_ACID, ρᵣ_ACID, ∂ρ∂pₗ_ACID, ∂ρ∂pᵣ_ACID, ∂Hₜ∂pₗ_ACID, ∂Hₜ∂pᵣ_ACID,
        Hₜₗ_ACID, Hₜᵣ_ACID, ∂ρ∂Tₗ_ACID, ∂ρ∂Tᵣ_ACID, ∂Hₜ∂Tₗ_ACID, ∂Hₜ∂Tᵣ_ACID,
        ∂ρ∂Y₁ₗ_ACID, ∂ρ∂Y₁ᵣ_ACID =
        EOS_ACID(
            pₗ,pᵣ,
            cells[face.owner].var[👉.u],cells[face.neighbour].var[👉.u],
            cells[face.owner].var[👉.v],cells[face.neighbour].var[👉.v],
            cells[face.owner].var[👉.w],cells[face.neighbour].var[👉.w],
            cells[face.owner].var[👉.T],cells[face.neighbour].var[👉.T],
            cells[face.owner].var[👉.Y₁],cells[face.neighbour].var[👉.Y₁]
            
        )

        #Y₁ₗ_ACID = Y₁ᵣ
        #Y₁ᵣ_ACID = Y₁ₗ

        ∂ρ∂Y₁ₗ_ACID = ∂ρ∂Y₁ₗ
        ∂ρ∂Y₁ᵣ_ACID = ∂ρ∂Y₁ᵣ
        Y₁ₗ_ACID = Y₁ₗ
        Y₁ᵣ_ACID = Y₁ᵣ
     
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( Wₗ * ρₗ * Uₙ * ΔS + Wₗ * ∂ρ∂Y₁ₗ * Y₁ₗ * Uₙ * ΔS )
        #push!(A_vals, ( Wᵣ * ρᵣ * Uₙ * ΔS + Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ᵣ * Uₙ * ΔS ))
        push!(A_vals, ( Wᵣ * ρᵣ_ACID * Uₙ * ΔS + Wᵣ * ∂ρ∂Y₁ᵣ_ACID * Y₁ᵣ_ACID * Uₙ * ΔS ))
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( Wᵣ * ρᵣ * Uₙ * ΔS + Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ᵣ * Uₙ * ΔS )
        #push!(A_vals, -( Wₗ * ρₗ * Uₙ * ΔS + Wₗ * ∂ρ∂Y₁ₗ * Y₁ₗ * Uₙ * ΔS ))
        push!(A_vals, -( Wₗ * ρₗ_ACID * Uₙ * ΔS + Wₗ * ∂ρ∂Y₁ₗ_ACID * Y₁ₗ_ACID * Uₙ * ΔS ))


        # convective terms
        #B[face.owner] -= (Wₗ * ρₗ * Y₁ₗ + Wᵣ * ρᵣ * Y₁ᵣ) * Uₙ * ΔS
        #B[face.neighbour] += (Wₗ * ρₗ * Y₁ₗ + Wᵣ * ρᵣ * Y₁ᵣ) * Uₙ * ΔS
        ρₙYₙ = ( Wₗ * ρₗ * 0.5 * ( Y₁ₗ + cells[face.owner].var[👉.Y₁ⁿ] ) + 
        Wᵣ * ρᵣ_ACID * 0.5 * ( Y₁ᵣ_ACID + cells[face.neighbour].var[👉.Y₁ⁿ] ) )
        B[face.owner] -= ρₙYₙ * Uₙ * ΔS
        ρₙYₙ = ( Wₗ * ρₗ_ACID * 0.5 * ( Y₁ₗ_ACID + cells[face.owner].var[👉.Y₁ⁿ] ) + 
        Wᵣ * ρᵣ * 0.5 * ( Y₁ᵣ + cells[face.neighbour].var[👉.Y₁ⁿ] ) )
        B[face.neighbour] += ρₙYₙ * Uₙ * ΔS
        
=#

    end
    

    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    bc_wall = []

    bc_slipwall = []
    append!( bc_slipwall, faces_boundary_top )
    append!( bc_slipwall, faces_boundary_bottom )
    
    bc_subinlet = []
    #append!( bc_subinlet, faces_boundary_left )
    
    bc_suboutlet = []
    #append!( bc_supoutlet, faces_boundary_right )
    
    bc_supoutlet = []
    append!( bc_supoutlet, faces_boundary_left )
    append!( bc_supoutlet, faces_boundary_right )
    
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
