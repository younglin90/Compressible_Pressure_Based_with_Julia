

function momentum!(
    š::controls,
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
        A_vals[diagon] = cell.var[š.Ļ]*cell.Ī©/š.Īt

        B[diagon, 1] = -( cell.var[š.Ļ]*cell.var[š.u] - cell.var[š.Ļāæ]*cell.var[š.uāæ] )*cell.Ī©/š.Īt
            
        B[diagon, 2] = -( cell.var[š.Ļ]*cell.var[š.v] - cell.var[š.Ļāæ]*cell.var[š.vāæ] )*cell.Ī©/š.Īt

        # gravity
        B[diagon, 2] += cell.var[š.Ļ]*cell.Ī© * (-9.8)

        diagon += 1

    end

    
    āĪpāx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pā = 0.5 * (cells[face.owner].var[š.p] + cells[face.neighbour].var[š.p])
        āĪpāx0[face.owner, 1] += pā * face.nĢ[1] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 2] += pā * face.nĢ[2] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 3] += pā * face.nĢ[3] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.neighbour, 1] -= pā * face.nĢ[1] * face.ĪS / cells[face.neighbour].Ī©
        āĪpāx0[face.neighbour, 2] -= pā * face.nĢ[2] * face.ĪS / cells[face.neighbour].Ī©
        āĪpāx0[face.neighbour, 3] -= pā * face.nĢ[3] * face.ĪS / cells[face.neighbour].Ī©
    end

    for face in faces_boundary
        pā = cells[face.owner].var[š.p]
        āĪpāx0[face.owner, 1] += pā * face.nĢ[1] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 2] += pā * face.nĢ[2] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 3] += pā * face.nĢ[3] * face.ĪS / cells[face.owner].Ī©
    end

    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal

        Ļā = cells[face.owner].var[š.Ļ]
        Ļįµ£ = cells[face.neighbour].var[š.Ļ]
        pā = cells[face.owner].var[š.p]
        pįµ£ = cells[face.neighbour].var[š.p]
        uā = cells[face.owner].var[š.u]
        uįµ£ = cells[face.neighbour].var[š.u]
        vā = cells[face.owner].var[š.v]
        vįµ£ = cells[face.neighbour].var[š.v]
        Ī¼ā = cells[face.owner].var[š.Ī¼]
        Ī¼įµ£ = cells[face.neighbour].var[š.Ī¼]
        Uāā = uā * face.nĢ[1] + vā * face.nĢ[2]
        Uāįµ£ = uįµ£ * face.nĢ[1] + vįµ£ * face.nĢ[2]
        Uā = 0.5 * (Uāā + Uāįµ£)
        ĪS = face.ĪS

        centerā = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerįµ£ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ĪLR = norm(centerįµ£ - centerā)

        invĻĪt = 0.5 * (1.0/Ļā + 1.0/Ļįµ£) * š.Īt
        
        # Rhie-Chow
        Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 1] * face.nĢ[1]
        Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 2] * face.nĢ[2]
        #Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 3] * face.nĢ[3]
        Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 1] * face.nĢ[1]
        Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 2] * face.nĢ[2]
        #Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 3] * face.nĢ[3]
        Uā -= invĻĪt * (pįµ£-pā) / ĪLR

        Wā = 0.5 * (1.0 + sign(Uā))
        Wįµ£ = 1.0 - Wā

        Ļā = Wā * Ļā + Wįµ£ * Ļįµ£
        uā = Wā * uā + Wįµ£ * uįµ£
        vā = Wā * vā + Wįµ£ * vįµ£
        Ī¼ā = Wā * Ī¼ā + Wįµ£ * Ī¼įµ£


        #--- ACID ----
        Ļā_ACID, Ļįµ£_ACID, āĻāpā_ACID, āĻāpįµ£_ACID, āHāāpā_ACID, āHāāpįµ£_ACID,
        Hāā_ACID, Hāįµ£_ACID, āĻāTā_ACID, āĻāTįµ£_ACID, āHāāTā_ACID, āHāāTįµ£_ACID,
        āĻāYāā_ACID, āĻāYāįµ£_ACID =
        EOS_ACID_vf(
            pā,pįµ£,
            cells[face.owner].var[š.u],cells[face.neighbour].var[š.u],
            cells[face.owner].var[š.v],cells[face.neighbour].var[š.v],
            cells[face.owner].var[š.w],cells[face.neighbour].var[š.w],
            cells[face.owner].var[š.T],cells[face.neighbour].var[š.T],
            cells[face.owner].var[š.Ī±ā],cells[face.neighbour].var[š.Ī±ā]
            
        )

        #Ļā_ACID = Ļā
        #Ļįµ£_ACID = Ļįµ£

        Ļā = Wā * Ļā + Wįµ£ * Ļįµ£
        
        #Ļā = ( Wā * Ļā + Wįµ£ * Ļįµ£_ACID )
        #Ļā = Ļā
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( Ļā * Wā * Uā * ĪS )
        push!(A_vals, ( Ļā * Wįµ£ * Uā * ĪS ))
        
        #Ļā = ( Wā * Ļā_ACID + Wįµ£ * Ļįµ£ )
        #Ļā = Ļįµ£
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( Ļā * Wįµ£ * Uā * ĪS )
        push!(A_vals, -( Ļā * Wā * Uā * ĪS ))


        # convective terms
        #Ļā = ( Wā * Ļā + Wįµ£ * Ļįµ£_ACID )
        #Ļā = Ļā
        B[face.owner, 1] -= Ļā * uā * Uā * ĪS
        B[face.neighbour, 1] += Ļā * uā * Uā * ĪS

        #Ļā = ( Wā * Ļā_ACID + Wįµ£ * Ļįµ£ )
        #Ļā = Ļįµ£
        B[face.owner, 2] -= Ļā * vā * Uā * ĪS
        B[face.neighbour, 2] += Ļā * vā * Uā * ĪS

        # pressure terms
        pā = 0.5 * (pā + pįµ£)

        B[face.owner, 1] -= pā * face.nĢ[1] * ĪS
        B[face.neighbour, 1] += pā * face.nĢ[1] * ĪS 
        
        B[face.owner, 2] -= pā * face.nĢ[2] * ĪS
        B[face.neighbour, 2] += pā * face.nĢ[2] * ĪS

#=
        # viscous terms
        B[face.owner, 1] += Ī¼ā * (uįµ£ - uā) / ĪLR * ĪS
        B[face.neighbour, 1] -= Ī¼ā * (uįµ£ - uā) / ĪLR * ĪS
        
        B[face.owner, 2] += Ī¼ā * (vįµ£ - vā) / ĪLR * ĪS
        B[face.neighbour, 2] -= Ī¼ā * (vįµ£ - vā) / ĪLR * ĪS
=#

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        Ļā = cells[face.owner].var[š.Ļ]
        ĪS = face.ĪS

        Uā = 0.0
        Uā += cells[face.owner].var[š.u]*face.nĢ[1]
        Uā += cells[face.owner].var[š.v]*face.nĢ[2]
        Uā += cells[face.owner].var[š.w]*face.nĢ[3]
        Uā0 = Uā

        invU = cells[face.owner].var[š.u] - Uā * face.nĢ[1]
        invV = cells[face.owner].var[š.v] - Uā * face.nĢ[2]
        invW = cells[face.owner].var[š.w] - Uā * face.nĢ[3]

        Uā = invU * face.nĢ[1]
        Uā += invV * face.nĢ[2]
        Uā += invW * face.nĢ[3]
        
        #A_vals[face.owner] += Ļā * Uā * ĪS

        Uā = 0.0

        # convective terms
        B[face.owner, 1] -= Ļā * invU * Uā * ĪS
        B[face.owner, 2] -= Ļā * invV * Uā * ĪS
        
        # pressure terms
        pā = cells[face.owner].var[š.p]
        B[face.owner, 1] -= pā * face.nĢ[1] * ĪS
        B[face.owner, 2] -= pā * face.nĢ[2] * ĪS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ĪU = solve(ps, A, B)
    




    relax = 0.7




    diagon = 1
    maximum_U = -1.e12
    for cell in cells

        cell.var[š.u] += relax*ĪU[diagon, 1]
        cell.var[š.v] += relax*ĪU[diagon, 2]
        
        maximum_U = max(maximum_U,abs(cell.var[š.u]))
        maximum_U = max(maximum_U,abs(cell.var[š.v]))
        maximum_U = max(maximum_U,abs(cell.var[š.w]))

        diagon += 1
    end


    #return log10(norm(ĪU))
    return log10(norm(ĪU)/(maximum_U+1.e-20))
   

end



#=
function momentum!(
    š::controls,
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
        A_vals[diagon] = cell.var[š.Ļ]*cell.Ī©/š.Īt

        B[diagon, 1] = -( cell.var[š.Ļ]*cell.var[š.u] - cell.var[š.Ļāæ]*cell.var[š.uāæ] )*cell.Ī©/š.Īt
            
        B[diagon, 2] = -( cell.var[š.Ļ]*cell.var[š.v] - cell.var[š.Ļāæ]*cell.var[š.vāæ] )*cell.Ī©/š.Īt

        # gravity
        B[diagon, 2] += cell.var[š.Ļ]*cell.Ī© * (-9.8)

        diagon += 1

    end

    
    āĪpāx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pā = 0.5 * (cells[face.owner].var[š.p] + cells[face.neighbour].var[š.p])
        āĪpāx0[face.owner, 1] += pā * face.nĢ[1] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 2] += pā * face.nĢ[2] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 3] += pā * face.nĢ[3] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.neighbour, 1] -= pā * face.nĢ[1] * face.ĪS / cells[face.neighbour].Ī©
        āĪpāx0[face.neighbour, 2] -= pā * face.nĢ[2] * face.ĪS / cells[face.neighbour].Ī©
        āĪpāx0[face.neighbour, 3] -= pā * face.nĢ[3] * face.ĪS / cells[face.neighbour].Ī©
    end

    for face in faces_boundary
        pā = cells[face.owner].var[š.p]
        āĪpāx0[face.owner, 1] += pā * face.nĢ[1] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 2] += pā * face.nĢ[2] * face.ĪS / cells[face.owner].Ī©
        āĪpāx0[face.owner, 3] += pā * face.nĢ[3] * face.ĪS / cells[face.owner].Ī©
    end

    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal

        Ļā = cells[face.owner].var[š.Ļ]
        Ļįµ£ = cells[face.neighbour].var[š.Ļ]
        pā = cells[face.owner].var[š.p]
        pįµ£ = cells[face.neighbour].var[š.p]
        uā = cells[face.owner].var[š.u]
        uįµ£ = cells[face.neighbour].var[š.u]
        vā = cells[face.owner].var[š.v]
        vįµ£ = cells[face.neighbour].var[š.v]
        Ī¼ā = cells[face.owner].var[š.Ī¼]
        Ī¼įµ£ = cells[face.neighbour].var[š.Ī¼]
        Uāā = uā * face.nĢ[1] + vā * face.nĢ[2]
        Uāįµ£ = uįµ£ * face.nĢ[1] + vįµ£ * face.nĢ[2]
        Uā = 0.5 * (Uāā + Uāįµ£)
        ĪS = face.ĪS

        centerā = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerįµ£ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ĪLR = norm(centerįµ£ - centerā)

        invĻĪt = 0.5 * (1.0/Ļā + 1.0/Ļįµ£) * š.Īt
        
        # Rhie-Chow
        Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 1] * face.nĢ[1]
        Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 2] * face.nĢ[2]
        Uā += 0.5 * š.Īt / Ļā * āĪpāx0[face.owner, 3] * face.nĢ[3]
        Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 1] * face.nĢ[1]
        Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 2] * face.nĢ[2]
        Uā += 0.5 * š.Īt / Ļįµ£ * āĪpāx0[face.neighbour, 3] * face.nĢ[3]
        Uā -= invĻĪt * (pįµ£-pā) / ĪLR

        Wā = 0.5 * (1.0 + sign(Uā))
        Wįµ£ = 1.0 - Wā

        Ļā = Wā * Ļā + Wįµ£ * Ļįµ£
        uā = Wā * uā + Wįµ£ * uįµ£
        vā = Wā * vā + Wįµ£ * vįµ£
        Ī¼ā = Wā * Ī¼ā + Wįµ£ * Ī¼įµ£


        #--- ACID ----
        Ļā_ACID, Ļįµ£_ACID, āĻāpā_ACID, āĻāpįµ£_ACID, āHāāpā_ACID, āHāāpįµ£_ACID,
        Hāā_ACID, Hāįµ£_ACID, āĻāTā_ACID, āĻāTįµ£_ACID, āHāāTā_ACID, āHāāTįµ£_ACID,
        āĻāYāā_ACID, āĻāYāįµ£_ACID =
        EOS_ACID_vf(
            pā,pįµ£,
            cells[face.owner].var[š.u],cells[face.neighbour].var[š.u],
            cells[face.owner].var[š.v],cells[face.neighbour].var[š.v],
            cells[face.owner].var[š.w],cells[face.neighbour].var[š.w],
            cells[face.owner].var[š.T],cells[face.neighbour].var[š.T],
            cells[face.owner].var[š.Ī±ā],cells[face.neighbour].var[š.Ī±ā]
            
        )

        #Ļā_ACID = Ļā
        #Ļįµ£_ACID = Ļįµ£

        Ļā = Wā * Ļā + Wįµ£ * Ļįµ£
        
        Ļā = ( Wā * Ļā + Wįµ£ * Ļįµ£_ACID )
        #Ļā = Ļā
        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)
        A_vals[face.owner] += ( Ļā * Wā * Uā * ĪS )
        push!(A_vals, ( Ļā * Wįµ£ * Uā * ĪS ))
        
        Ļā = ( Wā * Ļā_ACID + Wįµ£ * Ļįµ£ )
        #Ļā = Ļįµ£
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)
        A_vals[face.neighbour] -= ( Ļā * Wįµ£ * Uā * ĪS )
        push!(A_vals, -( Ļā * Wā * Uā * ĪS ))


        # convective terms
        Ļā = ( Wā * Ļā + Wįµ£ * Ļįµ£_ACID )
        #Ļā = Ļā
        B[face.owner, 1] -= Ļā * uā * Uā * ĪS
        B[face.neighbour, 1] += Ļā * uā * Uā * ĪS

        Ļā = ( Wā * Ļā_ACID + Wįµ£ * Ļįµ£ )
        #Ļā = Ļįµ£
        B[face.owner, 2] -= Ļā * vā * Uā * ĪS
        B[face.neighbour, 2] += Ļā * vā * Uā * ĪS

        # pressure terms
        pā = 0.5 * (pā + pįµ£)

        B[face.owner, 1] -= pā * face.nĢ[1] * ĪS
        B[face.neighbour, 1] += pā * face.nĢ[1] * ĪS 
        
        B[face.owner, 2] -= pā * face.nĢ[2] * ĪS
        B[face.neighbour, 2] += pā * face.nĢ[2] * ĪS

#=
        # viscous terms
        B[face.owner, 1] += Ī¼ā * (uįµ£ - uā) / ĪLR * ĪS
        B[face.neighbour, 1] -= Ī¼ā * (uįµ£ - uā) / ĪLR * ĪS
        
        B[face.owner, 2] += Ī¼ā * (vįµ£ - vā) / ĪLR * ĪS
        B[face.neighbour, 2] -= Ī¼ā * (vįµ£ - vā) / ĪLR * ĪS
=#

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        Ļā = cells[face.owner].var[š.Ļ]
        ĪS = face.ĪS

        Uā = 0.0
        Uā += cells[face.owner].var[š.u]*face.nĢ[1]
        Uā += cells[face.owner].var[š.v]*face.nĢ[2]
        Uā += cells[face.owner].var[š.w]*face.nĢ[3]
        Uā0 = Uā

        invU = cells[face.owner].var[š.u] - Uā * face.nĢ[1]
        invV = cells[face.owner].var[š.v] - Uā * face.nĢ[2]
        invW = cells[face.owner].var[š.w] - Uā * face.nĢ[3]

        Uā = invU * face.nĢ[1]
        Uā += invV * face.nĢ[2]
        Uā += invW * face.nĢ[3]
        
        #A_vals[face.owner] += Ļā * Uā * ĪS

        Uā = 0.0

        # convective terms
        B[face.owner, 1] -= Ļā * invU * Uā * ĪS
        B[face.owner, 2] -= Ļā * invV * Uā * ĪS
        
        # pressure terms
        pā = cells[face.owner].var[š.p]
        B[face.owner, 1] -= pā * face.nĢ[1] * ĪS
        B[face.owner, 2] -= pā * face.nĢ[2] * ĪS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ĪU = solve(ps, A, B)
    




    relax = 0.1




    diagon = 1
    maximum_U = -1.e12
    for cell in cells

        cell.var[š.u] += relax*ĪU[diagon, 1]
        cell.var[š.v] += relax*ĪU[diagon, 2]
        
        maximum_U = max(maximum_U,abs(cell.var[š.u]))
        maximum_U = max(maximum_U,abs(cell.var[š.v]))
        maximum_U = max(maximum_U,abs(cell.var[š.w]))

        diagon += 1
    end


    #return log10(norm(ĪU))
    return log10(norm(ĪU)/(maximum_U+1.e-20))
   

end
=#