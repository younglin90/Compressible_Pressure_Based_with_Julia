

function momentum!(
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

    # contruct A matrix diagonal terms
    # contruct B vector
    A_rows = zeros(Int64, length(cells))
    A_cols = zeros(Int64, length(cells))
    A_vals = zeros(Float64, length(cells))
    B = zeros(Float64, length(cells), 2)
    
    diagon = 1

    for cell in cells
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon
        A_vals[diagon] = cell.var[👉.ρ]*cell.Ω/👉.Δt

        B[diagon, 1] = -( cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ] )*cell.Ω/👉.Δt
            
        B[diagon, 2] = -( cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ] )*cell.Ω/👉.Δt

        # gravity
        B[diagon, 2] += cell.var[👉.ρ]*cell.Ω * (-9.8)

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
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        invρΔt = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        #Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        #Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        μₙ = Wₗ * μₗ + Wᵣ * μᵣ


        #--- ACID ----
        ρₗ_ACID, ρᵣ_ACID, ∂ρ∂pₗ_ACID, ∂ρ∂pᵣ_ACID, ∂Hₜ∂pₗ_ACID, ∂Hₜ∂pᵣ_ACID,
        Hₜₗ_ACID, Hₜᵣ_ACID, ∂ρ∂Tₗ_ACID, ∂ρ∂Tᵣ_ACID, ∂Hₜ∂Tₗ_ACID, ∂Hₜ∂Tᵣ_ACID,
        ∂ρ∂Y₁ₗ_ACID, ∂ρ∂Y₁ᵣ_ACID =
        EOS_ACID_vf(
            pₗ,pᵣ,
            cells[face.owner].var[👉.u],cells[face.neighbour].var[👉.u],
            cells[face.owner].var[👉.v],cells[face.neighbour].var[👉.v],
            cells[face.owner].var[👉.w],cells[face.neighbour].var[👉.w],
            cells[face.owner].var[👉.T],cells[face.neighbour].var[👉.T],
            cells[face.owner].var[👉.α₁],cells[face.neighbour].var[👉.α₁]
            
        )

        #ρₗ_ACID = ρₗ
        #ρᵣ_ACID = ρᵣ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        #ρₙ = ( Wₗ * ρₗ + Wᵣ * ρᵣ_ACID )
        #ρₙ = ρₗ
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( ρₙ * Wₗ * Uₙ * ΔS )
        push!(A_vals, ( ρₙ * Wᵣ * Uₙ * ΔS ))
        
        #ρₙ = ( Wₗ * ρₗ_ACID + Wᵣ * ρᵣ )
        #ρₙ = ρᵣ
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( ρₙ * Wᵣ * Uₙ * ΔS )
        push!(A_vals, -( ρₙ * Wₗ * Uₙ * ΔS ))


        # convective terms
        #ρₙ = ( Wₗ * ρₗ + Wᵣ * ρᵣ_ACID )
        #ρₙ = ρₗ
        B[face.owner, 1] -= ρₙ * uₙ * Uₙ * ΔS
        B[face.neighbour, 1] += ρₙ * uₙ * Uₙ * ΔS

        #ρₙ = ( Wₗ * ρₗ_ACID + Wᵣ * ρᵣ )
        #ρₙ = ρᵣ
        B[face.owner, 2] -= ρₙ * vₙ * Uₙ * ΔS
        B[face.neighbour, 2] += ρₙ * vₙ * Uₙ * ΔS

        # pressure terms
        pₙ = 0.5 * (pₗ + pᵣ)

        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.neighbour, 1] += pₙ * face.n̂[1] * ΔS 
        
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        B[face.neighbour, 2] += pₙ * face.n̂[2] * ΔS

#=
        # viscous terms
        B[face.owner, 1] += μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        B[face.neighbour, 1] -= μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        
        B[face.owner, 2] += μₙ * (vᵣ - vₗ) / ΔLR * ΔS
        B[face.neighbour, 2] -= μₙ * (vᵣ - vₗ) / ΔLR * ΔS
=#

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ρₗ = cells[face.owner].var[👉.ρ]
        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]
        Uₙ0 = Uₙ

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]
        
        #A_vals[face.owner] += ρₗ * Uₙ * ΔS

        Uₙ = 0.0

        # convective terms
        B[face.owner, 1] -= ρₗ * invU * Uₙ * ΔS
        B[face.owner, 2] -= ρₗ * invV * Uₙ * ΔS
        
        # pressure terms
        pₙ = cells[face.owner].var[👉.p]
        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ΔU = solve(ps, A, B)
    




    relax = 0.7




    diagon = 1
    maximum_U = -1.e12
    for cell in cells

        cell.var[👉.u] += relax*ΔU[diagon, 1]
        cell.var[👉.v] += relax*ΔU[diagon, 2]
        
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))
        maximum_U = max(maximum_U,abs(cell.var[👉.w]))

        diagon += 1
    end


    #return log10(norm(ΔU))
    return log10(norm(ΔU)/(maximum_U+1.e-20))
   

end



#=
function momentum!(
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

    # contruct A matrix diagonal terms
    # contruct B vector
    A_rows = zeros(Int64, length(cells))
    A_cols = zeros(Int64, length(cells))
    A_vals = zeros(Float64, length(cells))
    B = zeros(Float64, length(cells), 2)
    
    diagon = 1

    for cell in cells
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon
        A_vals[diagon] = cell.var[👉.ρ]*cell.Ω/👉.Δt

        B[diagon, 1] = -( cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ] )*cell.Ω/👉.Δt
            
        B[diagon, 2] = -( cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ] )*cell.Ω/👉.Δt

        # gravity
        B[diagon, 2] += cell.var[👉.ρ]*cell.Ω * (-9.8)

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
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        invρΔt = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        μₙ = Wₗ * μₗ + Wᵣ * μᵣ


        #--- ACID ----
        ρₗ_ACID, ρᵣ_ACID, ∂ρ∂pₗ_ACID, ∂ρ∂pᵣ_ACID, ∂Hₜ∂pₗ_ACID, ∂Hₜ∂pᵣ_ACID,
        Hₜₗ_ACID, Hₜᵣ_ACID, ∂ρ∂Tₗ_ACID, ∂ρ∂Tᵣ_ACID, ∂Hₜ∂Tₗ_ACID, ∂Hₜ∂Tᵣ_ACID,
        ∂ρ∂Y₁ₗ_ACID, ∂ρ∂Y₁ᵣ_ACID =
        EOS_ACID_vf(
            pₗ,pᵣ,
            cells[face.owner].var[👉.u],cells[face.neighbour].var[👉.u],
            cells[face.owner].var[👉.v],cells[face.neighbour].var[👉.v],
            cells[face.owner].var[👉.w],cells[face.neighbour].var[👉.w],
            cells[face.owner].var[👉.T],cells[face.neighbour].var[👉.T],
            cells[face.owner].var[👉.α₁],cells[face.neighbour].var[👉.α₁]
            
        )

        #ρₗ_ACID = ρₗ
        #ρᵣ_ACID = ρᵣ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        ρₙ = ( Wₗ * ρₗ + Wᵣ * ρᵣ_ACID )
        #ρₙ = ρₗ
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( ρₙ * Wₗ * Uₙ * ΔS )
        push!(A_vals, ( ρₙ * Wᵣ * Uₙ * ΔS ))
        
        ρₙ = ( Wₗ * ρₗ_ACID + Wᵣ * ρᵣ )
        #ρₙ = ρᵣ
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( ρₙ * Wᵣ * Uₙ * ΔS )
        push!(A_vals, -( ρₙ * Wₗ * Uₙ * ΔS ))


        # convective terms
        ρₙ = ( Wₗ * ρₗ + Wᵣ * ρᵣ_ACID )
        #ρₙ = ρₗ
        B[face.owner, 1] -= ρₙ * uₙ * Uₙ * ΔS
        B[face.neighbour, 1] += ρₙ * uₙ * Uₙ * ΔS

        ρₙ = ( Wₗ * ρₗ_ACID + Wᵣ * ρᵣ )
        #ρₙ = ρᵣ
        B[face.owner, 2] -= ρₙ * vₙ * Uₙ * ΔS
        B[face.neighbour, 2] += ρₙ * vₙ * Uₙ * ΔS

        # pressure terms
        pₙ = 0.5 * (pₗ + pᵣ)

        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.neighbour, 1] += pₙ * face.n̂[1] * ΔS 
        
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        B[face.neighbour, 2] += pₙ * face.n̂[2] * ΔS

#=
        # viscous terms
        B[face.owner, 1] += μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        B[face.neighbour, 1] -= μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        
        B[face.owner, 2] += μₙ * (vᵣ - vₗ) / ΔLR * ΔS
        B[face.neighbour, 2] -= μₙ * (vᵣ - vₗ) / ΔLR * ΔS
=#

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ρₗ = cells[face.owner].var[👉.ρ]
        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]
        Uₙ0 = Uₙ

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]
        
        #A_vals[face.owner] += ρₗ * Uₙ * ΔS

        Uₙ = 0.0

        # convective terms
        B[face.owner, 1] -= ρₗ * invU * Uₙ * ΔS
        B[face.owner, 2] -= ρₗ * invV * Uₙ * ΔS
        
        # pressure terms
        pₙ = cells[face.owner].var[👉.p]
        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ΔU = solve(ps, A, B)
    




    relax = 0.1




    diagon = 1
    maximum_U = -1.e12
    for cell in cells

        cell.var[👉.u] += relax*ΔU[diagon, 1]
        cell.var[👉.v] += relax*ΔU[diagon, 2]
        
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))
        maximum_U = max(maximum_U,abs(cell.var[👉.w]))

        diagon += 1
    end


    #return log10(norm(ΔU))
    return log10(norm(ΔU)/(maximum_U+1.e-20))
   

end
=#