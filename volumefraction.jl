

function volumefraction_advection!(
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

    A_rows::Vector{Int64} = []
    A_cols::Vector{Int64} = []
    A_vals::Vector{Float64} = []
    
    # contruct A matrix diagonal terms
    # contruct B vector
    B = zeros(Float64, length(cells))
    
    diagon = 1

    for cell in cells
        
        push!(A_rows, diagon)
        push!(A_cols, diagon)

        tmp_A_var = cell.Ω/👉.Δt

        push!(A_vals, tmp_A_var)

        B[diagon] = -(cell.var[👉.α₁]-cell.var[👉.α₁ⁿ])*cell.Ω/👉.Δt

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
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        #=
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        =#
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR

        wₗ = 0.5 * (1.0 + sign(Uₙ))
        wᵣ = 1.0 - wₗ
        
        uₙ = wₗ * uₗ + wᵣ * uᵣ
        vₙ = wₗ * vₗ + wᵣ * vᵣ
        
        α₁ₗ = cells[face.owner].var[👉.α₁]
        α₁ᵣ = cells[face.neighbour].var[👉.α₁]
        α₁ₙ = wₗ * α₁ₗ + wᵣ * α₁ᵣ

        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += wₗ * Uₙ * ΔS
        push!(A_vals, wᵣ * Uₙ * ΔS)
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= wᵣ * Uₙ * ΔS
        push!(A_vals, -wₗ * Uₙ * ΔS)

        # convective terms
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        B[face.neighbour] += α₁ₙ * Uₙ * ΔS

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
        α₁ₙ = cells[face.owner].var[👉.α₁]
        
        A_vals[face.owner] += Uₙ * ΔS

        # convective terms
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        
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

        α₁ₙ = cells[face.owner].var[👉.α₁]
        
        A_vals[face.owner] += Uₙ * ΔS

        # convective terms
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        
    end
    
    for face in bc_subinlet
        ΔS = face.ΔS

        uₙ = 1.0
        vₙ = 0.0

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        α₁ₙ = 1.0

        #A_vals[face.owner] += Uₙ * ΔS

        # convective terms
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        
    end

    for face in bc_suboutlet
        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        α₁ₙ = cells[face.owner].var[👉.α₁]

        A_vals[face.owner] += Uₙ * ΔS

        # convective terms
        α₁ₙ = cells[face.owner].var[👉.α₁]
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        
    end
    

    for face in bc_supoutlet
        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]

        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        α₁ₙ = cells[face.owner].var[👉.α₁]

        A_vals[face.owner] += Uₙ * ΔS

        # convective terms
        α₁ₙ = cells[face.owner].var[👉.α₁]
        B[face.owner] -= α₁ₙ * Uₙ * ΔS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    Δα₁ = solve(ps, A, B)
    
    #Δα₁ = gmres(A, B)
    



    relax = 0.9



    diagon = 1
    for cell in cells

        cell.var[👉.α₁] += relax*Δα₁[diagon]

        diagon += 1
    end


    return norm(Δα₁)
   

end


#=

function volumefraction_advection!(
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
        
        A_rows[diagon] = diagon
        A_cols[diagon] = diagon
        A_vals[diagon] = cell.Ω/👉.Δt

        B[diagon] = -(cell.var[👉.α₁] - cell.var[👉.α₁ⁿ])*cell.Ω/👉.Δt


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

    UₙA = zeros(Float64, length(cells))

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

        invρΔt = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        #=
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR
        =#

        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ
        
        α₁ₗ = cells[face.owner].var[👉.α₁]
        α₁ᵣ = cells[face.neighbour].var[👉.α₁]
        α₁ₙ = Wₗ * α₁ₗ + Wᵣ * α₁ᵣ
     
        α₁ₗⁿ = cells[face.owner].var[👉.α₁ⁿ]
        α₁ᵣⁿ = cells[face.neighbour].var[👉.α₁ⁿ]
        α₁ₙⁿ = Wₗ * α₁ₗⁿ + Wᵣ * α₁ᵣⁿ

#=
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += 0.5 * Wₗ * Uₙ * ΔS
        push!(A_vals, 0.5 * Wᵣ * Uₙ * ΔS)
        #A_vals[face.owner] += 0.5 * Uₙₗ * ΔS
        #push!(A_vals, 0.5 * Uₙₗ * ΔS)
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= 0.5 * Wᵣ * Uₙ * ΔS
        push!(A_vals, -0.5 * Wₗ * Uₙ * ΔS)
        #A_vals[face.neighbour] -= 0.5 * Uₙᵣ * ΔS
        #push!(A_vals, -0.5 * Uₙᵣ * ΔS)
=#

        # convective terms
        #B[face.owner] -= 0.5 * ( α₁ₙ + α₁ₙⁿ ) * Uₙ * ΔS
        #B[face.neighbour] += 0.5 * ( α₁ₙ + α₁ₙⁿ ) * Uₙ * ΔS
        B[face.owner] -= (uₗ * face.n̂[1] + vₗ * face.n̂[2]) * ( Wₗ * α₁ₗ + Wᵣ * α₁ᵣ ) * ΔS
        B[face.neighbour] += (uᵣ * face.n̂[1] + vᵣ * face.n̂[2]) * ( Wₗ * α₁ₗ + Wᵣ * α₁ᵣ ) * ΔS
        

#=
        # conv
        Kpₗ = 0.0
        Kpᵣ = 0.0
        B[face.owner] += ( 0.5*( α₁ₗ + α₁ₗⁿ ) + Kpₗ ) * Uₙ * ΔS
        B[face.neighbour] -= ( 0.5*( α₁ᵣ + α₁ᵣⁿ ) + Kpᵣ ) * Uₙ * ΔS

        #UₙA[face.owner] += Uₙ * ΔS
        #UₙA[face.neighbour] -= Uₙ * ΔS
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

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]
        
        Uₙ = 0.0

        #A_vals[face.owner] += ρₗ * Uₙ * ΔS

        # convective terms
        #Y₁ₙ = cells[face.owner].var[👉.Y₁]
        #B[face.owner] -= Y₁ₙ * ρₗ * Uₙ * ΔS
        
    end
 



    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    Δα₁ = solve(ps, A, B)
    



    relax = 0.1


    diagon = 1
    for cell in cells

        cell.var[👉.α₁] += relax*Δα₁[diagon]

        diagon += 1
    end


    return log10(norm(Δα₁))
   

end
=#


function volumefraction!(
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

    A_rows::Vector{Int64} = []
    A_cols::Vector{Int64} = []
    A_vals::Vector{Float64} = []
    
    # contruct A matrix diagonal terms
    # contruct B vector
    B = zeros(Float64, length(cells))
    
    diagon = 1

    for cell in cells
        
        push!(A_rows, diagon)
        push!(A_cols, diagon)

        push!(A_vals, cell.Ω/👉.Δt)

        B[diagon] = -(
            cell.var[👉.α₁] - cell.var[👉.α₁ⁿ]
            )*cell.Ω/👉.Δt

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
        
        α₁ₗ = cells[face.owner].var[👉.α₁]
        α₁ᵣ = cells[face.neighbour].var[👉.α₁]
        α₁ₙ = Wₗ * α₁ₗ + Wᵣ * α₁ᵣ
     
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        push!(A_vals, Wᵣ * Uₙ * ΔS)
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        push!(A_vals, -Wₗ * Uₙ * ΔS)

        A_vals[face.owner] += Wₗ * Uₙ * ΔS
        A_vals[face.neighbour] -= Wᵣ * Uₙ * ΔS


        # convective terms
        flux = (Wₗ * α₁ₗ + Wᵣ * α₁ᵣ) * Uₙ * ΔS
        B[face.owner] -= flux
        B[face.neighbour] += flux

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

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]
        
        Uₙ = 0.0

        A_vals[face.owner] += ρₗ * Uₙ * ΔS

        # convective terms
        Y₁ₙ = cells[face.owner].var[👉.Y₁]
        B[face.owner] -= Y₁ₙ * ρₗ * Uₙ * ΔS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    Δα₁ = solve(ps, A, B)
    





    relax = 0.9




    diagon = 1
    for cell in cells

        cell.var[👉.α₁] += relax*Δα₁[diagon]

        diagon += 1
    end


    return log10(norm(Δα₁))
   

end
