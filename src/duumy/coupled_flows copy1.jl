
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

        Ω = cell.Ω
        Δt = 👉.Δt
        u = cell.var[👉.u]
        v = cell.var[👉.v]
        ρ = cell.var[👉.ρ]
        Hₜ = cell.var[👉.Hₜ]
        p = cell.var[👉.p]
        ∂Hₜ∂p = cell.var[👉.∂Hₜ∂p]
        ∂Hₜ∂T = cell.var[👉.∂Hₜ∂T]
        ∂ρ∂p = cell.var[👉.∂ρ∂p]
        ∂ρ∂T = cell.var[👉.∂ρ∂T]
        ρⁿ = cell.var[👉.ρⁿ]
        uⁿ = cell.var[👉.uⁿ]
        vⁿ = cell.var[👉.vⁿ]
        Hₜⁿ = cell.var[👉.Hₜⁿ]
        pⁿ = cell.var[👉.pⁿ]

        #println(pⁿ,uⁿ,vⁿ,Hₜⁿ)
        
        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        A_vals[i] = ∂ρ∂p*Ω/Δt
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        i += 1
        A_rows[i] = ijStart + 1;  A_cols[i] = ijStart + 4
        A_vals[i] = ∂ρ∂T*Ω/Δt

        B[ijStart + 1] = -(ρ - ρⁿ)*Ω/Δt

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1
        A_vals[i] = ∂ρ∂p*Ω/Δt * u

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = ρ*Ω/Δt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 4
        A_vals[i] = ∂ρ∂T*Ω/Δt * u

        B[ijStart + 2] = -(ρ*u - ρⁿ*uⁿ)*cell.Ω/Δt

        # y-momentum
        #g = -9.8
        g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        A_vals[i] = ∂ρ∂p*Ω/Δt * v + ∂ρ∂p*Ω * g
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = ρ*Ω/Δt

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 4
        A_vals[i] = ∂ρ∂T*Ω/Δt * v + ∂ρ∂T*Ω * g

        B[ijStart + 3] = -(ρ*v - ρⁿ*vⁿ)*Ω/Δt + ρ*Ω * g



        # energy
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 1
        A_vals[i] = ∂ρ∂p*Ω/Δt * Hₜ + ∂Hₜ∂p*Ω/Δt * ρ - Ω/Δt

        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 2
        A_vals[i] = u*Ω/Δt * ρ
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 3
        A_vals[i] = v*Ω/Δt * ρ
        
        i += 1
        A_rows[i] = ijStart + 4; A_cols[i] = ijStart + 4
        A_vals[i] = ∂ρ∂T*Ω/Δt * Hₜ + ∂Hₜ∂T*Ω/Δt * ρ
        
        B[ijStart + 4] = -(ρ*Hₜ - ρⁿ*Hₜⁿ)*Ω/Δt + (p - pⁿ)*Ω/Δt


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
        #Y₁ₗ = cells[face.owner].var[👉.Y₁]
        #Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]
        #∂ρ∂Y₁ₗ = cells[face.owner].var[👉.∂ρ∂Y₁]
        #∂ρ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂ρ∂Y₁]
        #∂Hₜ∂Y₁ₗ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        #∂Hₜ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂Hₜ∂Y₁]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d = 0.5 * (cells[face.owner].Ω / (Ap[face.owner]+1.e-250) + cells[face.neighbour].Ω / (Ap[face.neighbour]+1.e-250) )
        #d̂ = 👉.Δt / ρˢ
        d̂ = d / (1.0 + ρˢ / 👉.Δt * d)
        
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

        ∂ρ∂pₗ_ACID = ∂ρ∂pₗ
        ∂ρ∂pᵣ_ACID = ∂ρ∂pᵣ
        ∂Hₜ∂pₗ_ACID = ∂Hₜ∂pₗ
        ∂Hₜ∂pᵣ_ACID = ∂Hₜ∂pᵣ
        ∂ρ∂Tₗ_ACID = ∂ρ∂Tₗ
        ∂ρ∂Tᵣ_ACID = ∂ρ∂Tᵣ
        ∂Hₜ∂Tₗ_ACID = ∂Hₜ∂Tₗ
        ∂Hₜ∂Tᵣ_ACID = ∂Hₜ∂Tᵣ

        ρₗ_ACID = ρₗ
        ρᵣ_ACID = ρᵣ
        Hₜₗ_ACID = Hₜₗ
        Hₜᵣ_ACID = Hₜᵣ

        #=
        ρₙₗ_ACID = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        ρₙᵣ_ACID = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        
        Hₜₙₗ_ACID = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        Hₜₙᵣ_ACID = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        =#

#=


        #------------------------
        # continuity
        
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙₗ_ACID * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Uₙ * ΔS - ρₙₗ_ACID * d̂ / ΔLR * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙᵣ_ACID * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Uₙ * ΔS + ρₙᵣ_ACID * d̂ / ΔLR * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[1] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[2] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Uₙ * ΔS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙₗ_ACID * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙₗ_ACID * d̂ / ΔLR * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙᵣ_ACID * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙᵣ_ACID * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[2] * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * uₙ * Uₙ * ΔS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] +=  ( Wₗ * ∂ρ∂pₗ       * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙₗ_ACID * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙₗ_ACID * d̂ / ΔLR * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ      * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙᵣ_ACID * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙᵣ_ACID * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ_ACID * vₙ * Uₙ *ΔS ))


        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hₜₙₗ_ACID * Uₙ * ΔS + ρₙₗ_ACID * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙₗ_ACID * d̂ / ΔLR * Hₜₙₗ_ACID * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Hₜₙₗ_ACID * Uₙ * ΔS + ρₙₗ_ACID * Wᵣ * ∂Hₜ∂pᵣ_ACID * Uₙ * ΔS - ρₙₗ_ACID * d̂ / ΔLR * Hₜₙₗ_ACID * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hₜₙᵣ_ACID * Uₙ * ΔS + ρₙᵣ_ACID * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙᵣ_ACID * d̂ / ΔLR * Hₜₙᵣ_ACID * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Hₜₙᵣ_ACID * Uₙ * ΔS + ρₙᵣ_ACID * Wₗ * ∂Hₜ∂pₗ_ACID * Uₙ * ΔS + ρₙᵣ_ACID * d̂ / ΔLR * Hₜₙᵣ_ACID * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[1] * Hₜₙₗ_ACID * ΔS + ρₙₗ_ACID * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[1] * Hₜₙₗ_ACID * ΔS + ρₙₗ_ACID * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[1] * Hₜₙᵣ_ACID * ΔS + ρₙᵣ_ACID * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[1] * Hₜₙᵣ_ACID * ΔS + ρₙᵣ_ACID * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙₗ_ACID * 0.5 * face.n̂[2] * Hₜₙₗ_ACID * ΔS + ρₙₗ_ACID * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙₗ_ACID * 0.5 * face.n̂[2] * Hₜₙₗ_ACID * ΔS + ρₙₗ_ACID * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙᵣ_ACID * 0.5 * face.n̂[2] * Hₜₙᵣ_ACID * ΔS + ρₙᵣ_ACID * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙᵣ_ACID * 0.5 * face.n̂[2] * Hₜₙᵣ_ACID * ΔS + ρₙᵣ_ACID * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙₗ_ACID * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Hₜᵣ_ACID * Uₙ * ΔS + ρₙₗ_ACID * Wᵣ * ∂Hₜ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ_ACID * Hₜᵣ * Uₙ * ΔS + ρₙᵣ_ACID * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Hₜₗ_ACID * Uₙ * ΔS + ρₙᵣ_ACID * Wₗ * ∂Hₜ∂Tₗ_ACID * Uₙ * ΔS ))
        

        # B
        B[ijStartₗ + 1] -= ρₙₗ_ACID * Uₙ * ΔS
        B[ijStartᵣ + 1] += ρₙᵣ_ACID * Uₙ * ΔS
        # B
        B[ijStartₗ + 2] -= ( ρₙₗ_ACID * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( ρₙᵣ_ACID * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        # B
        B[ijStartₗ + 3] -= ( ρₙₗ_ACID * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( ρₙᵣ_ACID * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        # B
        B[ijStartₗ + 4] -= ρₙₗ_ACID * Hₜₙₗ_ACID * Uₙ * ΔS
        B[ijStartᵣ + 4] += ρₙᵣ_ACID * Hₜₙᵣ_ACID * Uₙ * ΔS

        


=#


        #------------------------
        # continuity
        # p'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        #ρᵣ = ρᵣ_ACID
        #∂ρ∂pᵣ = ∂ρ∂pᵣ_ACID
        #∂Hₜ∂pᵣ = ∂Hₜ∂pᵣ_ACID
        #Hₜᵣ = Hₜᵣ_ACID
        #∂ρ∂Tᵣ = ∂ρ∂Tᵣ_ACID
        #∂Hₜ∂Tᵣ = ∂Hₜ∂Tᵣ_ACID

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Uₙ * ΔS - ρₙ * d̂ / ΔLR * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        #ρₗ = ρₗ_ACID
        #∂ρ∂pₗ = ∂ρ∂pₗ_ACID
        #∂Hₜ∂pₗ = ∂Hₜ∂pₗ_ACID
        #Hₜₗ = Hₜₗ_ACID
        #∂ρ∂Tₗ = ∂ρ∂Tₗ_ACID
        #∂Hₜ∂Tₗ = ∂Hₜ∂Tₗ_ACID

        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Uₙ * ΔS + ρₙ * d̂ / ΔLR * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        
        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Uₙ * ΔS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * d̂ / ΔLR * uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ_ACID * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ_ACID * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * uₙ * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * uₙ * Uₙ * ΔS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] +=  ( Wₗ * ∂ρ∂pₗ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * d̂ / ΔLR * vₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ_ACID * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ_ACID * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * vₙ * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ_ACID * vₙ * Uₙ *ΔS ))


        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * Hₜₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ_ACID * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ_ACID * Uₙ * ΔS - ρₙ * d̂ / ΔLR * Hₜₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * Hₜₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ_ACID * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ_ACID * Uₙ * ΔS + ρₙ * d̂ / ΔLR * Hₜₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ_ACID * Hₜᵣ_ACID * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ_ACID * Uₙ * ΔS ))
        
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ_ACID * Hₜₗ_ACID * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ_ACID * Uₙ * ΔS ))
        


        # B

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        B[ijStartᵣ + 1] += ( ρₙ * Uₙ * ΔS )
        # B

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        B[ijStartᵣ + 2] += ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        # B

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        B[ijStartᵣ + 3] += ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        # B

        ρₙ = (Wₗ*ρₗ + Wᵣ*ρᵣ_ACID)
        Hₜₙ = (Wₗ * Hₜₗ + Wᵣ * Hₜᵣ_ACID)
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        ρₙ = (Wₗ*ρₗ_ACID + Wᵣ*ρᵣ)
        Hₜₙ = (Wₗ * Hₜₗ_ACID + Wᵣ * Hₜᵣ)
        B[ijStartᵣ + 4] += ( ρₙ * Hₜₙ * Uₙ * ΔS )




#=

        #------------------------
        # continuity
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * ΔS ))

        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * d̂ / ΔLR * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS ))


        

        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] +=  ( Wₗ * ∂ρ∂pₗ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * d̂ / ΔLR * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS ))


        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * Hₜₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * Hₜₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * d̂ / ΔLR * Hₜₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * d̂ / ΔLR * Hₜₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS ))


        # B
        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartᵣ + 1] += ( ρₙ * Uₙ * ΔS )
        # B
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        # B
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        # B
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartᵣ + 4] += ( ρₙ * Hₜₙ * Uₙ * ΔS )
=#


    end


    
    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )

    bc_slipwall = []
    append!( bc_slipwall, faces_boundary_top )
    append!( bc_slipwall, faces_boundary_bottom )
    
    bc_subinlet = []
    append!( bc_subinlet, faces_boundary_left )
    
    bc_suboutlet = []
    append!( bc_suboutlet, faces_boundary_right )

    for face in bc_slipwall
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]

        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        uₙ = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        vₙ = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        wₙ = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]


        ρₙ, Hₜₙ, cₙ = faceEOS_vf!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        # continuity
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS

        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        

    end
 
    for face in bc_subinlet
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]

        ΔS = face.ΔS

        uₙ = 1.0
        #uₙ = 0.5 * ( 1.0 + cells[face.owner].var[👉.u] )
        vₙ = 0.0
        wₙ = 0.0
        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        Tₙ = 300.0
        #Tₙ = 0.5 * ( 300.0 + cells[face.owner].var[👉.T] )
        α₁ₙ = 1.0
        
        #pₙ = 101325.0

        ρₙ, Hₜₙ, cₙ = faceEOS_vf!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = 1.0 * norm(centerᵣ - centerₗ)

        # continuity
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * ΔS# + ρₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * ΔS)

        
        # x-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * uₙ * Uₙ * ΔS# + ρₙ * uₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * uₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * uₙ * Uₙ * ΔS)

        
        # y-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * vₙ * Uₙ * ΔS# + ρₙ * vₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * ρₙ * vₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * vₙ * Uₙ * ΔS)


        # energy
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS)


        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        

    end
 
    for face in bc_suboutlet
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]

        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]
        wₙ = cells[face.owner].var[👉.w]
        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]

        pₙ = 0.5 * ( 101325.0 + cells[face.owner].var[👉.p] )
        
        ρₙ, Hₜₙ, cₙ = faceEOS_vf!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = 2.0 * norm(centerᵣ - centerₗ)

        # continuity
        i += 1
        A_vals[i] += 0.5 * (∂ρ∂pₙ * Uₙ * ΔS) + ρₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += 0.5 * (∂ρ∂pₙ * uₙ * Uₙ * ΔS) + ρₙ * uₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += ρₙ * uₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += 0.5 * (∂ρ∂pₙ * vₙ * Uₙ * ΔS) + ρₙ * vₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * vₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += 0.5 * (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS) + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS

        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        

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

    relax_p = 1.0
    relax_U = 1.0
    relax_T = 1.0


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

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]

        cell.var[👉.p] = max(cell.var[👉.p],0.1)
        cell.var[👉.T] = max(cell.var[👉.T],0.1)
        
        #println(cell.var[👉.p])
        
        norm_p += ΔQ[ijStart + 1]^2
        norm_U += ΔQ[ijStart + 2]^2
        norm_U += ΔQ[ijStart + 3]^2
        norm_T += ΔQ[ijStart + 4]^2
        maximum_p = max(maximum_p,abs(cell.var[👉.p]))
        maximum_U = max(maximum_U,abs(cell.var[👉.u]))
        maximum_U = max(maximum_U,abs(cell.var[👉.v]))
        maximum_T = max(maximum_T,abs(cell.var[👉.T]))
        maximum_Y = max(maximum_Y,abs(cell.var[👉.Y₁]))

        diagon += 1
    end

    #sleep(1000.0)


    return norm(ΔQ),log10(sqrt(norm_p)/length(cells)/(maximum_p+1.e-200)),
    log10(sqrt(norm_U)/length(cells)/(maximum_U+1.e-200)),
    log10(sqrt(norm_T)/length(cells)/(maximum_T+1.e-200))
   

end
