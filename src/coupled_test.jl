
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

    # contruct A matrix diagonal terms
    # contruct B vector
   # N = length(cells) * 5

    #A_rows::Vector{Int64} = []
    #A_cols::Vector{Int64} = []
    #A_vals::Vector{Float64} = []

    #B::Vector{Float64} = []
    
    B_n = 4
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
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt

        B[ijStart + 1] =
        -(cell.var[👉.ρ] - cell.var[👉.ρⁿ])*cell.Ω/👉.Δt

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.u]
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.u]

        B[ijStart + 2] = 
        -(cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ])*cell.Ω/👉.Δt

        # y-momentum
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = 
        cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.v] +
        cell.var[👉.∂ρ∂p]*cell.Ω * (-9.8)
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = 
        cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.v] +
        cell.var[👉.∂ρ∂T]*cell.Ω * (-9.8)

        B[ijStart + 3] = 
        -(cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ])*cell.Ω/👉.Δt +
        cell.var[👉.ρ]*cell.Ω * (-9.8)



        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = 
        cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.Hₜ] +
        cell.var[👉.∂Hₜ∂p]*cell.Ω/👉.Δt * cell.var[👉.ρ] -
        cell.Ω/👉.Δt

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = cell.var[👉.u]*cell.Ω/👉.Δt * cell.var[👉.ρ]
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = cell.var[👉.v]*cell.Ω/👉.Δt * cell.var[👉.ρ]
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = 
        cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.Hₜ] +
        cell.var[👉.∂Hₜ∂T]*cell.Ω/👉.Δt * cell.var[👉.ρ]

        B[ijStart + 4] = 
        -(cell.var[👉.ρ]*cell.var[👉.Hₜ] - cell.var[👉.ρⁿ]*cell.var[👉.Hₜⁿ])*cell.Ω/👉.Δt +
        (cell.var[👉.p] - cell.var[👉.pⁿ])*cell.Ω/👉.Δt



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
        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR

        flux = ρₙ * Uₙ * face.ΔS
        Ap[face.owner] += flux
        Ap[face.neighbour] -= flux
    end

    for face in faces_boundary
        #pₙ = cells[face.owner].var[👉.p]
        #Ap[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end



    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        wₗ = 0.0#cells[face.owner].var[👉.w]
        wᵣ = 0.0#cells[face.neighbour].var[👉.w]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        ∂ρ∂pₗ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂pᵣ = cells[face.neighbour].var[👉.∂ρ∂p]
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = cells[face.neighbour].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d = 0.5 * (cells[face.owner].Ω / (Ap[face.owner]+1.e-50) 
        + cells[face.neighbour].Ω / (Ap[face.neighbour]+1.e-50) )
        d̂ = 👉.Δt / ρˢ
        #d̂ = d / (1.0 + ρˢ / 👉.Δt * d)
        #d̂ = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR
        #=
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ -= d̂ * 👉.Δt * (pᵣ-pₗ) / ΔLR
        =#

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Hₜₙ = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ

        pₙ = 0.5 * (pₗ + pᵣ)

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)


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

        
        ∂ρ∂pₗ_ACID = ∂ρ∂pₗ
        ∂ρ∂pᵣ_ACID = ∂ρ∂pᵣ
        ∂Hₜ∂pₗ_ACID = ∂Hₜ∂pₗ
        ∂Hₜ∂pᵣ_ACID = ∂Hₜ∂pᵣ
        ∂ρ∂Tₗ_ACID = ∂ρ∂Tₗ
        ∂ρ∂Tᵣ_ACID = ∂ρ∂Tᵣ
        ∂Hₜ∂Tₗ_ACID = ∂Hₜ∂Tₗ
        ∂Hₜ∂Tᵣ_ACID = ∂Hₜ∂Tᵣ
        #=
        ρₗ_ACID = ρₗ
        ρᵣ_ACID = ρᵣ
        Hₜₗ_ACID = Hₜₗ
        Hₜᵣ_ACID = Hₜᵣ
        =#

        #------------------------
        # continuity
        
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * Uₙ * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Uₙ * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * ΔS
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * ΔS
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * ΔS
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * ΔS
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Uₙ * ΔS ))


        # B
        B[ijStartₗ + 1] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * Uₙ * ΔS
        B[ijStartᵣ + 1] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * Uₙ * ΔS
        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * uₙ * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * uₙ * ΔS
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * uₙ * ΔS
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * uₙ * Uₙ * ΔS ))


        # B
        B[ijStartₗ + 2] -= ( (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        

        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * vₙ * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * vₙ * ΔS
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * vₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * vₙ * ΔS
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ_ACID * vₙ * Uₙ *ΔS ))


        # B
        B[ijStartₗ + 3] -= ( (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hm = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hm * Uₙ * ΔS +
        ρm * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS +
        ρm * d̂ / ΔLR * Hm * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Hm * Uₙ * ΔS +
        ρm * Wᵣ * ∂Hₜ∂pᵣ_ACID * Uₙ * ΔS -
        ρm * d̂ / ΔLR * Hm * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hm * Uₙ * ΔS +
        ρm * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS -
        ρm * d̂ / ΔLR * Hm * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Hm * Uₙ * ΔS +
        ρm * Wₗ * ∂Hₜ∂pₗ_ACID * Uₙ * ΔS +
        ρm * d̂ / ΔLR * Hm * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hm = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID
        A_vals[iₗ] += ( ρm * 0.5 * face.n̂[1] * Hm * ΔS + ρm * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρm * 0.5 * face.n̂[1] * Hm * ΔS + ρm * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ
        A_vals[iᵣ] -= ( ρm * 0.5 * face.n̂[1] * Hm * ΔS + ρm * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρm * 0.5 * face.n̂[1] * Hm * ΔS + ρm * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hm = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID
        A_vals[iₗ] += ( ρm * 0.5 * face.n̂[2] * Hm * ΔS + ρm * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρm * 0.5 * face.n̂[2] * Hm * ΔS + ρm * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ
        A_vals[iᵣ] -= ( ρm * 0.5 * face.n̂[2] * Hm * ΔS + ρm * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρm * 0.5 * face.n̂[2] * Hm * ΔS + ρm * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hm = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρm * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ_ACID * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ
        A_vals[iᵣ] -= Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Hₜₗ_ACID * Uₙ * ΔS + ρm * Wₗ * ∂Hₜ∂Tₗ_ACID * Uₙ * ΔS ))

        # B
        B[ijStartₗ + 4] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID) * Uₙ * ΔS
        B[ijStartᵣ + 4] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ) * Uₙ * ΔS





    end
    
    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        pₙ = cells[face.owner].var[👉.p]
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
        

        Uₙ = 0.0

        # continuity
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * ΔS

        B[ijStartₗ + 1] -= 0.0
        
        # x-momentum
        i += 1
        A_vals[i] += face.n̂[1] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.u] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.u] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.u] * ΔS

        B[ijStartₗ + 2] -= pₙ * face.n̂[1] * ΔS
        
        # y-momentum
        i += 1
        A_vals[i] += face.n̂[2] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS

        B[ijStartₗ + 3] -= pₙ * face.n̂[2] * ΔS

        # energy
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS

        B[ijStartₗ + 4] -= 0.0
        
    end
 


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ΔQ = solve(ps, A, B)
    #ΔQ = A\B

   # error()
   # exit()

   #println(length(B)," ",length(ΔQ))
   #println(norm(A*ΔQ-B))

    relax_p = 0.3
    relax_U = 0.5
    relax_T = 0.3

    diagon = 1
    maximum_p = -1.e12
    maximum_U = -1.e12
    maximum_T = -1.e12
    norm_p = 0.0
    norm_U = 0.0
    norm_T = 0.0
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]
        
        #println(cell.var[👉.p])
        
        norm_p += ΔQ[ijStart + 1]^2
        norm_U += ΔQ[ijStart + 2]^2
        norm_U += ΔQ[ijStart + 3]^2
        norm_T += ΔQ[ijStart + 4]^2
        maximum_p = max(maximum_p,abs(cell.var[👉.p]))
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))
        maximum_T = max(maximum_T,abs(cell.var[👉.T]))

        diagon += 1
    end

    #sleep(1000.0)


    return log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-20)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-20)),
    log10(sqrt(norm_T)/length(cells)/(maximum_T+1.e-20))
   

end



#=

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

    # contruct A matrix diagonal terms
    # contruct B vector
   # N = length(cells) * 5

    #A_rows::Vector{Int64} = []
    #A_cols::Vector{Int64} = []
    #A_vals::Vector{Float64} = []

    #B::Vector{Float64} = []
    
    B_n = 3
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
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3

        B[ijStart + 1] =
        -(cell.var[👉.ρ] - cell.var[👉.ρⁿ])*cell.Ω/👉.Δt

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.u]
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3
        
        B[ijStart + 2] = 
        -(cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ])*cell.Ω/👉.Δt

        # y-momentum
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.v]
        A_vals[i] += cell.var[👉.∂ρ∂p]*cell.Ω * (-9.8)
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt

        B[ijStart + 3] = 
        -(cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ])*cell.Ω/👉.Δt
        B[ijStart + 3] += cell.var[👉.ρ]*cell.Ω * (-9.8)



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
        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR

        flux = ρₙ * Uₙ * face.ΔS
        Ap[face.owner] += flux
        Ap[face.neighbour] -= flux
    end

    for face in faces_boundary
        #pₙ = cells[face.owner].var[👉.p]
        #Ap[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end



    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        wₗ = 0.0#cells[face.owner].var[👉.w]
        wᵣ = 0.0#cells[face.neighbour].var[👉.w]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        ∂ρ∂pₗ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂pᵣ = cells[face.neighbour].var[👉.∂ρ∂p]
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = cells[face.neighbour].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d = 0.5 * (cells[face.owner].Ω / (Ap[face.owner]+1.e-50) 
        + cells[face.neighbour].Ω / (Ap[face.neighbour]+1.e-50) )
        d̂ = 👉.Δt / ρˢ
        #d̂ = d / (1.0 + ρˢ / 👉.Δt * d)
        #d̂ = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR
        #=
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ -= d̂ * 👉.Δt * (pᵣ-pₗ) / ΔLR
        =#

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Hₜₙ = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ

        pₙ = 0.5 * (pₗ + pᵣ)

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)


        #--- ACID ----
        ρₗ_ACID,ρᵣ_ACID,∂ρ∂pₗ_ACID,∂ρ∂pᵣ_ACID =
        EOS_ACID(
            pₗ,pᵣ,
            cells[face.owner].var[👉.u],cells[face.neighbour].var[👉.u],
            cells[face.owner].var[👉.v],cells[face.neighbour].var[👉.v],
            cells[face.owner].var[👉.w],cells[face.neighbour].var[👉.w],
            cells[face.owner].var[👉.T],cells[face.neighbour].var[👉.T],
            cells[face.owner].var[👉.Y₁],cells[face.neighbour].var[👉.Y₁]
            
        )

        ∂ρ∂pₗ_ACID = ∂ρ∂pₗ
        ∂ρ∂pᵣ_ACID = ∂ρ∂pᵣ

        #------------------------
        # continuity
        
        # p'
        #convfluxₗ = Wₗ * ∂ρ∂pₗ * Uₙ * ΔS
        #convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS
        #difffluxₗ = ρₙ * d̂ / ΔLR * ΔS
        #difffluxᵣ = -ρₙ * d̂ / ΔLR * ΔS
        
        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 1, ijStartᵣ + 1,
        #iᵣ, ijStartᵣ + 1, ijStartₗ + 1,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * Uₙ * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Uₙ * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * ΔS ))
        
        # u'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[1] * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 1, ijStartᵣ + 2,
        #iᵣ, ijStartᵣ + 1, ijStartₗ + 2,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * ΔS
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * ΔS
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[2] * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 1, ijStartᵣ + 3,
        #iᵣ, ijStartᵣ + 1, ijStartₗ + 3,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * ΔS
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * ΔS
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * ΔS ))
        
        # B
        B[ijStartₗ + 1] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * Uₙ * ΔS
        B[ijStartᵣ + 1] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * Uₙ * ΔS
        

        #------------------------
        # x-momentum

        # p'
        #convfluxₗ = Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS
        #convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS
        #difffluxₗ = ρₙ * d̂ / ΔLR * uₙ * ΔS
        #difffluxᵣ = -ρₙ * d̂ / ΔLR * uₙ * ΔS

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 2, ijStartᵣ + 1,
        #iᵣ, ijStartᵣ + 2, ijStartₗ + 1,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * uₙ * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₙ * Uₙ * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρₙ * Uₙ * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 2, ijStartᵣ + 2,
        #iᵣ, ijStartᵣ + 2, ijStartₗ + 2,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        # v'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 2, ijStartᵣ + 3,
        #iᵣ, ijStartᵣ + 2, ijStartₗ + 3,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * uₙ * ΔS
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * uₙ * ΔS
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # B
        B[ijStartₗ + 2] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS
        B[ijStartᵣ + 2] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS
        

        #------------------------
        # y-momentum
        
        # p'
        #convfluxₗ = Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS
        #convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS
        #difffluxₗ = ρₙ * d̂ / ΔLR * vₙ * ΔS
        #difffluxᵣ = -ρₙ * d̂ / ΔLR * vₙ * ΔS

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 3, ijStartᵣ + 1,
        #iᵣ, ijStartᵣ + 3, ijStartₗ + 1,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, Wᵣ * ∂ρ∂pᵣ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * d̂ / ΔLR * vₙ * ΔS)
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 3, ijStartᵣ + 2,
        #iᵣ, ijStartᵣ + 3, ijStartₗ + 2,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * vₙ * ΔS
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[1] * vₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * vₙ * ΔS
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        #convfluxₗ = ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₙ * Uₙ * ΔS
        #convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρₙ * Uₙ * ΔS
        #difffluxₗ = 0.0
        #difffluxᵣ = 0.0

        #iₗ += 1; iᵣ += 1
        #push_A_conv_diff!(A_rows, A_cols, A_vals, 
        #iₗ, ijStartₗ + 3, ijStartᵣ + 3,
        #iᵣ, ijStartᵣ + 3, ijStartₗ + 3,
        #convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS)
        
        A_vals[iᵣ] -= (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))

        # B
        B[ijStartₗ + 3] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS
        B[ijStartᵣ + 3] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS
        
       # println(ρₗ," ",ρᵣ_ACID)
       # println(ρₗ_ACID," ",ρᵣ)


    end
    
    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ijStartₗ = B_n*(face.owner-1)

        iₗ = A_n*(face.owner-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        pₙ = cells[face.owner].var[👉.p]
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
        

        Uₙ = 0.0

        # continuity
        iₗ += 1
        A_vals[iₗ] += 0.0
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * ΔS

        B[ijStartₗ + 1] -= 0.0
        
        # x-momentum
        iₗ += 1
        A_vals[iₗ] += face.n̂[1] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.u] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.u] * ΔS

        B[ijStartₗ + 2] -= pₙ * face.n̂[1] * ΔS
        
        # y-momentum
        iₗ += 1
        A_vals[iₗ] += face.n̂[2] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.v] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS

        B[ijStartₗ + 3] -= pₙ * face.n̂[2] * ΔS
        
    end
 


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ΔQ = solve(ps, A, B)
    #ΔQ = A\B

   # error()
   # exit()

   #println(length(B)," ",length(ΔQ))
   #println(norm(A*ΔQ-B))

    relax_p = 0.5
    relax_U = 0.9
    relax_T = 0.0

    diagon = 1
    maximum_p = -1.e12
    maximum_U = -1.e12
    norm_p = 0.0
    norm_U = 0.0
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        
        #println(cell.var[👉.p])
        
        norm_p += ΔQ[ijStart + 1]^2
        norm_U += ΔQ[ijStart + 2]^2
        norm_U += ΔQ[ijStart + 3]^2
        maximum_p = max(maximum_p,abs(cell.var[👉.p]))
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))

        diagon += 1
    end

    #sleep(1000.0)


    return log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-20)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-20)),
    0.0
   

end


=#


#=
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

    # contruct A matrix diagonal terms
    # contruct B vector
   # N = length(cells) * 5

    #A_rows::Vector{Int64} = []
    #A_cols::Vector{Int64} = []
    #A_vals::Vector{Float64} = []

    #B::Vector{Float64} = []
    
    B_n = 4
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
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt

        B[ijStart + 1] =
        -(cell.var[👉.ρ] - cell.var[👉.ρⁿ])*cell.Ω/👉.Δt

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.u]
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.u]
        
        B[ijStart + 2] = 
        -(cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ])*cell.Ω/👉.Δt

        # y-momentum
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.v]
        A_vals[i] += cell.var[👉.∂ρ∂p]*cell.Ω * (-9.8)
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = cell.var[👉.ρ]*cell.Ω/👉.Δt

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.v]
        A_vals[i] += cell.var[👉.∂ρ∂T]*cell.Ω * (-9.8)

        B[ijStart + 3] = 
        -(cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ])*cell.Ω/👉.Δt
        B[ijStart + 3] += cell.var[👉.ρ]*cell.Ω * (-9.8)


        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = 
        cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.Hₜ]
        + cell.var[👉.∂Hₜ∂p]*cell.Ω/👉.Δt * cell.var[👉.ρ]
        - cell.Ω/👉.Δt

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = cell.var[👉.u]*cell.Ω/👉.Δt * cell.var[👉.ρ]
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = cell.var[👉.v]*cell.Ω/👉.Δt * cell.var[👉.ρ]

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = 
        cell.var[👉.∂ρ∂T]*cell.Ω/👉.Δt * cell.var[👉.Hₜ]
        + cell.var[👉.∂Hₜ∂T]*cell.Ω/👉.Δt * cell.var[👉.ρ]

        B[ijStart + 4] = 
        -(cell.var[👉.ρ]*cell.var[👉.Hₜ] - cell.var[👉.ρⁿ]*cell.var[👉.Hₜⁿ])*cell.Ω/👉.Δt
        +(cell.var[👉.p] - cell.var[👉.pⁿ])*cell.Ω/👉.Δt



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
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d̂ = 👉.Δt / ρˢ
        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        #Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        #Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR

        flux = ρₙ * Uₙ * face.ΔS
        Ap[face.owner] += flux
        Ap[face.neighbour] -= flux
    end

    for face in faces_boundary
        #pₙ = cells[face.owner].var[👉.p]
        #Ap[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        #Ap[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end



    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal
        
        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        wₗ = 0.0#cells[face.owner].var[👉.w]
        wᵣ = 0.0#cells[face.neighbour].var[👉.w]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        ∂ρ∂pₗ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂pᵣ = cells[face.neighbour].var[👉.∂ρ∂p]
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = cells[face.neighbour].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d = 0.5 * (cells[face.owner].Ω / (Ap[face.owner]+1.e-50) 
        + cells[face.neighbour].Ω / (Ap[face.neighbour]+1.e-50) )
        d̂ = 👉.Δt / ρˢ
        #d̂ = d / (1.0 + ρˢ / 👉.Δt * d)
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Hₜₙ = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ

        pₙ = 0.5 * (pₗ + pᵣ)

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)

        #------------------------
        # continuity
        
        # p'
        convfluxₗ = Wₗ * ∂ρ∂pₗ * Uₙ * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS
        difffluxₗ = ρₙ * d̂ / ΔLR * ΔS
        difffluxᵣ = -ρₙ * d̂ / ΔLR * ΔS

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 1, ijStartᵣ + 1,
        iᵣ, ijStartᵣ + 1, ijStartₗ + 1,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # u'
        convfluxₗ = ρₙ * 0.5 * face.n̂[1] * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 1, ijStartᵣ + 2,
        iᵣ, ijStartᵣ + 1, ijStartₗ + 2,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # v'
        convfluxₗ = ρₙ * 0.5 * face.n̂[2] * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 1, ijStartᵣ + 3,
        iᵣ, ijStartᵣ + 1, ijStartₗ + 3,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # T'
        convfluxₗ = Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 1, ijStartᵣ + 4,
        iᵣ, ijStartᵣ + 1, ijStartₗ + 4,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # B
        covflux = ρₙ * Uₙ * ΔS
        B[ijStartₗ + 1] -= covflux
        B[ijStartᵣ + 1] += covflux
        

        #------------------------
        # x-momentum

        # p'
        convfluxₗ = Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS
        difffluxₗ = ρₙ * d̂ / ΔLR * uₙ * ΔS
        difffluxᵣ = -ρₙ * d̂ / ΔLR * uₙ * ΔS

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 2, ijStartᵣ + 1,
        iᵣ, ijStartᵣ + 2, ijStartₗ + 1,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # u'
        convfluxₗ = ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 2, ijStartᵣ + 2,
        iᵣ, ijStartᵣ + 2, ijStartₗ + 2,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # v'
        convfluxₗ = ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 2, ijStartᵣ + 3,
        iᵣ, ijStartᵣ + 2, ijStartₗ + 3,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        # T'
        convfluxₗ = Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 2, ijStartᵣ + 4,
        iᵣ, ijStartᵣ + 2, ijStartₗ + 4,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # B
        covflux = ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS
        B[ijStartₗ + 2] -= covflux
        B[ijStartᵣ + 2] += covflux
        

        #------------------------
        # y-momentum
        
        # p'
        convfluxₗ = Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS
        difffluxₗ = ρₙ * d̂ / ΔLR * vₙ * ΔS
        difffluxᵣ = -ρₙ * d̂ / ΔLR * vₙ * ΔS

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 3, ijStartᵣ + 1,
        iᵣ, ijStartᵣ + 3, ijStartₗ + 1,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # u'
        convfluxₗ = ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 3, ijStartᵣ + 2,
        iᵣ, ijStartᵣ + 3, ijStartₗ + 2,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # v'
        convfluxₗ = ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 3, ijStartᵣ + 3,
        iᵣ, ijStartᵣ + 3, ijStartₗ + 3,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        # T'
        convfluxₗ = Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 3, ijStartᵣ + 4,
        iᵣ, ijStartᵣ + 3, ijStartₗ + 4,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # B
        covflux = ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS
        B[ijStartₗ + 3] -= covflux
        B[ijStartᵣ + 3] += covflux
        

        #------------------------
        # energy
        convfluxₗ = Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS
        difffluxₗ = ρₙ * d̂ / ΔLR * Hₜₙ * ΔS
        difffluxᵣ = -ρₙ * d̂ / ΔLR * Hₜₙ * ΔS

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 4, ijStartᵣ + 1,
        iᵣ, ijStartᵣ + 4, ijStartₗ + 1,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # u'
        convfluxₗ = ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 4, ijStartᵣ + 2,
        iᵣ, ijStartᵣ + 4, ijStartₗ + 2,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # v'
        convfluxₗ = ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS
        convfluxᵣ = ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 4, ijStartᵣ + 3,
        iᵣ, ijStartᵣ + 4, ijStartₗ + 3,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)
        
        # T'
        convfluxₗ = Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + Wₗ * ρₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS
        convfluxᵣ = Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + Wᵣ * ρᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS
        #convfluxₗ = ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS
        #convfluxᵣ = ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS
        difffluxₗ = 0.0
        difffluxᵣ = 0.0

        iₗ += 1; iᵣ += 1
        push_A_conv_diff!(A_rows, A_cols, A_vals, 
        iₗ, ijStartₗ + 4, ijStartᵣ + 4,
        iᵣ, ijStartᵣ + 4, ijStartₗ + 4,
        convfluxₗ, difffluxₗ, convfluxᵣ, difffluxᵣ)

        # B
        covflux = ρₙ * Hₜₙ * Uₙ * ΔS
        B[ijStartₗ + 4] -= covflux
        B[ijStartᵣ + 4] += covflux



    end
    
    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ijStartₗ = B_n*(face.owner-1)

        iₗ = A_n*(face.owner-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        pₙ = cells[face.owner].var[👉.p]
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
        

        Uₙ = 0.0

        # continuity
        iₗ += 1
        A_vals[iₗ] += 0.0
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0

        B[ijStartₗ + 1] -= 0.0
        
        # x-momentum
        iₗ += 1
        A_vals[iₗ] += face.n̂[1] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.u] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.u] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0

        B[ijStartₗ + 2] -= pₙ * face.n̂[1] * ΔS
        
        # y-momentum
        iₗ += 1
        A_vals[iₗ] += face.n̂[2] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.v] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0

        B[ijStartₗ + 3] -= pₙ * face.n̂[2] * ΔS
        
        # energy
        iₗ += 1
        A_vals[iₗ] += 0.0
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[1] * Hₜₗ * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0#ρₗ * 0.5 * face.n̂[2] * Hₜₗ * ΔS
        iₗ += 1
        A_vals[iₗ] += 0.0

        B[ijStartₗ + 4] -= 0.0
        
    end
 


    A = sparse(A_rows,A_cols,A_vals)

    #spy(A, marker=".", markersize=1)
    #gui()
    #sleep(1000.0)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ΔQ = solve(ps, A, B)
    #ΔQ = A\B

   # error()
   # exit()

   #println(length(B)," ",length(ΔQ))
   #println(norm(A*ΔQ-B))

    relax_p = 0.01
    relax_U = 0.01
    relax_T = 0.01

    diagon = 1
    maximum_p = -1.e12
    maximum_U = -1.e12
    norm_p = 0.0
    norm_U = 0.0
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]
        
        #println(cell.var[👉.p])
        
        norm_p += ΔQ[ijStart + 1]^2
        norm_U += ΔQ[ijStart + 2]^2
        norm_U += ΔQ[ijStart + 3]^2
        maximum_p = max(maximum_p,abs(cell.var[👉.p]))
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))

        diagon += 1
    end

    #sleep(1000.0)


    return log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-20)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-20)),
    0.0
   

end
=#
