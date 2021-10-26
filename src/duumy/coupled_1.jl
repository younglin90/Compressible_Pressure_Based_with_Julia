
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
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 5
        A_vals[i] = cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt

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

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 5
        A_vals[i] = cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.u]

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

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 5
        A_vals[i] = 
        cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.v] +
        cell.var[👉.∂ρ∂Y₁]*cell.Ω * (-9.8)

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
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 5
        A_vals[i] = 
        cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.Hₜ] +
        cell.var[👉.∂Hₜ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.ρ]

        B[ijStart + 4] = 
        -(cell.var[👉.ρ]*cell.var[👉.Hₜ] - cell.var[👉.ρⁿ]*cell.var[👉.Hₜⁿ])*cell.Ω/👉.Δt +
        (cell.var[👉.p] - cell.var[👉.pⁿ])*cell.Ω/👉.Δt



        # mass fraction
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 1
        A_vals[i] = cell.var[👉.∂ρ∂p]*cell.Ω/👉.Δt * cell.var[👉.Y₁]

        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 4
        A_vals[i] = cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.Y₁]
        
        i += 1
        A_rows[i] = ijStart + 5; A_cols[i] = ijStart + 5
        A_vals[i] = 
        cell.var[👉.∂ρ∂Y₁]*cell.Ω/👉.Δt * cell.var[👉.Y₁] +
        cell.Ω/👉.Δt * cell.var[👉.ρ]

        B[ijStart + 5] = -(cell.var[👉.ρ]*cell.var[👉.Y₁] - cell.var[👉.ρⁿ]*cell.var[👉.Y₁ⁿ])*cell.Ω/👉.Δt




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
        Y₁ₗ = cells[face.owner].var[👉.Y₁]
        Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]
        ∂ρ∂Y₁ₗ = cells[face.owner].var[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ₗ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂Hₜ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂Hₜ∂Y₁]

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
        Y₁ₗ_ACID = Y₁ᵣ
        Y₁ᵣ_ACID = Y₁ₗ
        ∂Hₜ∂Y₁ₗ_ACID = ∂Hₜ∂Y₁ₗ
        ∂Hₜ∂Y₁ᵣ_ACID = ∂Hₜ∂Y₁ᵣ


        #=
        Y₁ₗ_ACID = Y₁ₗ
        Y₁ᵣ_ACID = Y₁ᵣ
        ∂ρ∂pₗ_ACID = ∂ρ∂pₗ
        ∂ρ∂pᵣ_ACID = ∂ρ∂pᵣ
        ∂Hₜ∂pₗ_ACID = ∂Hₜ∂pₗ
        ∂Hₜ∂pᵣ_ACID = ∂Hₜ∂pᵣ
        ∂ρ∂Tₗ_ACID = ∂ρ∂Tₗ
        ∂ρ∂Tᵣ_ACID = ∂ρ∂Tᵣ
        ∂Hₜ∂Tₗ_ACID = ∂Hₜ∂Tₗ
        ∂Hₜ∂Tᵣ_ACID = ∂Hₜ∂Tᵣ
        ∂ρ∂Y₁ₗ_ACID = ∂ρ∂Y₁ₗ
        ∂ρ∂Y₁ᵣ_ACID = ∂ρ∂Y₁ᵣ
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
        
        # Y'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ_ACID * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ_ACID * Uₙ * ΔS ))


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

        # Y'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ_ACID * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ_ACID * uₙ * Uₙ * ΔS ))

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

        # Y'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ_ACID * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -(  Wₗ * ∂ρ∂Y₁ₗ_ACID * vₙ * Uₙ *ΔS ))

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
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Hₜᵣ_ACID * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ
        A_vals[iᵣ] -= Wᵣ * ∂ρ∂Tᵣ_ACID * Hₜᵣ * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Hₜₗ_ACID * Uₙ * ΔS + ρm * Wₗ * ∂Hₜ∂Tₗ_ACID * Uₙ * ΔS ))

        
        # Y'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hm = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Hₜₗ * Uₙ * ΔS + ρm * Wₗ * ∂Hₜ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ_ACID * Hₜᵣ_ACID * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Y₁ᵣ_ACID * Uₙ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hm = Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ
        A_vals[iᵣ] -= Wᵣ * ∂ρ∂Y₁ᵣ * Hₜᵣ * Uₙ * ΔS + ρm * Wᵣ * ∂Hₜ∂Y₁ᵣ * Uₙ * ΔS
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ_ACID * Hₜₗ_ACID * Uₙ * ΔS + ρm * Wₗ * ∂Hₜ∂Y₁ₗ_ACID * Uₙ * ΔS ))

        # B
        B[ijStartₗ + 4] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID) * Uₙ * ΔS
        B[ijStartᵣ + 4] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ) * Uₙ * ΔS


        #------------------------
        # mass fraction

        #Y₁ₗ_ACID = Y₁ₗ
        #Y₁ᵣ_ACID = Y₁ᵣ

        # p'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Ym = (Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Ym * Uₙ * ΔS + ρm * d̂ / ΔLR * Ym * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Ym * Uₙ * ΔS + ρm * d̂ / ΔLR * Ym * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Ym = (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Ym * Uₙ * ΔS + ρm * d̂ / ΔLR * Ym * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Ym * Uₙ * ΔS + ρm * d̂ / ΔLR * Ym * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Ym = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID
        A_vals[iₗ] += ( ρm * 0.5 * face.n̂[1] * Ym * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρm * 0.5 * face.n̂[1] * Ym * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Ym = (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ)
        A_vals[iᵣ] -= ( ρm * 0.5 * face.n̂[1] * Ym * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρm * 0.5 * face.n̂[1] * Ym * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Ym = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID
        A_vals[iₗ] += ( ρm * 0.5 * face.n̂[2] * Ym * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρm * 0.5 * face.n̂[2] * Ym * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Ym = (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ)
        A_vals[iᵣ] -= ( ρm * 0.5 * face.n̂[2] * Ym * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρm * 0.5 * face.n̂[2] * Ym * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Ym = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Y₁ᵣ_ACID * Uₙ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Ym = (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Y₁ₗ_ACID * Uₙ * ΔS ))

        
        # Y'
        iₗ += 1; iᵣ += 1
        ρm = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Ym = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Y₁ₗ * Uₙ * ΔS + ρm * Wₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ_ACID * Y₁ᵣ_ACID * Uₙ * ΔS + ρm * Wᵣ * Uₙ * ΔS ))
        
        ρm = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Ym = (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ᵣ * Uₙ * ΔS + ρm * Wᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ_ACID * Y₁ₗ_ACID * Uₙ * ΔS + ρm * Wₗ * Uₙ * ΔS ))

        # B
        B[ijStartₗ + 5] -= (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID) * (Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ_ACID) * Uₙ * ΔS
        B[ijStartᵣ + 5] += (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ) * (Wₗ * Y₁ₗ_ACID + Wᵣ * Y₁ᵣ) * Uₙ * ΔS
        #B[ijStartₗ + 5] -= (Wₗ*ρₗ + Wᵣ*ρᵣ) * (Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ) * Uₙ * ΔS
        #B[ijStartᵣ + 5] += (Wₗ*ρₗ + Wᵣ*ρᵣ) * (Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ) * Uₙ * ΔS






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
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS

        B[ijStartₗ + 4] -= 0.0
        

        # mass fraction
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[1] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS
        i += 1
        A_vals[i] += 0.0#ρₗ * 0.5 * face.n̂[2] * cells[face.owner].var[👉.v] * ΔS

        B[ijStartₗ + 5] -= 0.0
        
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

    relax_p = 0.1
    relax_U = 0.3
    relax_T = 0.3
    relax_Y = 0.3

    diagon = 1
    maximum_p = -1.e12
    maximum_U = -1.e12
    maximum_T = -1.e12
    maximum_Y = -1.e12
    norm_p = 0.0
    norm_U = 0.0
    norm_T = 0.0
    norm_Y = 0.0
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]
        cell.var[👉.Y₁] += relax_Y * ΔQ[ijStart + 5]
        
        #println(cell.var[👉.p])
        
        norm_p += ΔQ[ijStart + 1]^2
        norm_U += ΔQ[ijStart + 2]^2
        norm_U += ΔQ[ijStart + 3]^2
        norm_T += ΔQ[ijStart + 4]^2
        norm_Y += ΔQ[ijStart + 5]^2
        maximum_p = max(maximum_p,abs(cell.var[👉.p]))
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))
        maximum_T = max(maximum_T,abs(cell.var[👉.T]))
        maximum_Y = max(maximum_Y,abs(cell.var[👉.Y₁]))

        diagon += 1
    end

    #sleep(1000.0)


    return log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-20)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-20)),
    log10(sqrt(norm_T)/length(cells)/(maximum_T+1.e-20))
   

end
