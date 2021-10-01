

function energy!(
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
    B = zeros(Float64, length(cells))
    
    diagon = 1

    for cell in cells
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon
        A_vals[diagon] = cell.var[👉.∂ρ∂T]*cell.var[👉.Hₜ]*cell.Ω/👉.Δt
        A_vals[diagon] += cell.var[👉.ρ]*cell.var[👉.∂Hₜ∂T]*cell.Ω/👉.Δt

        B[diagon] = -( 
            (cell.var[👉.ρ]*cell.var[👉.Hₜ] - cell.var[👉.p]) - 
            (cell.var[👉.ρⁿ]*cell.var[👉.Hₜⁿ] - cell.var[👉.pⁿ])
            )*cell.Ω/👉.Δt
        
        #println(A_vals[diagon]," " , B[diagon])

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
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
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
#=        
        tmp_A_var_convL = Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS
        tmp_A_var_convR = Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS
        
        tmp_A_var_convL += Wₗ * ρₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS
        tmp_A_var_convR += Wᵣ * ρᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS

        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        push!(A_vals, tmp_A_var_convR)
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        push!(A_vals, -tmp_A_var_convL)

        A_vals[face.owner] += tmp_A_var_convL
        A_vals[face.neighbour] -= tmp_A_var_convR
=#
        # convective terms
        covflux = ( Wₗ * ρₗ * Hₜₗ + Wᵣ * ρᵣ * Hₜᵣ ) * Uₙ * ΔS
        B[face.owner] -= covflux
        B[face.neighbour] += covflux

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

        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ΔT = solve(ps, A, B)
    




    relax = 0.5




    diagon = 1
    maximum_T = -1.e12
    for cell in cells

        cell.var[👉.T] += relax*ΔT[diagon]

        println(cell.var[👉.T])

        maximum_T = max(maximum_T,cell.var[👉.T])

        diagon += 1
    end


    return log10(norm(ΔT)/(maximum_T+1.e-20))
   

end
