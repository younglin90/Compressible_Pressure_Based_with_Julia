

function energy!(
    ð::controls,
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
        A_vals[diagon] = cell.var[ð.âÏâT]*cell.var[ð.Hâ]*cell.Î©/ð.Ît
        A_vals[diagon] += cell.var[ð.Ï]*cell.var[ð.âHââT]*cell.Î©/ð.Ît

        B[diagon] = -( 
            (cell.var[ð.Ï]*cell.var[ð.Hâ] - cell.var[ð.p]) - 
            (cell.var[ð.Ïâ¿]*cell.var[ð.Hââ¿] - cell.var[ð.pâ¿])
            )*cell.Î©/ð.Ît
        
        #println(A_vals[diagon]," " , B[diagon])

        diagon += 1

    end

    
    âÎpâx0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pâ = 0.5 * (cells[face.owner].var[ð.p] + cells[face.neighbour].var[ð.p])
        âÎpâx0[face.owner, 1] += pâ * face.nÌ[1] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 2] += pâ * face.nÌ[2] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 3] += pâ * face.nÌ[3] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.neighbour, 1] -= pâ * face.nÌ[1] * face.ÎS / cells[face.neighbour].Î©
        âÎpâx0[face.neighbour, 2] -= pâ * face.nÌ[2] * face.ÎS / cells[face.neighbour].Î©
        âÎpâx0[face.neighbour, 3] -= pâ * face.nÌ[3] * face.ÎS / cells[face.neighbour].Î©
    end

    for face in faces_boundary
        pâ = cells[face.owner].var[ð.p]
        âÎpâx0[face.owner, 1] += pâ * face.nÌ[1] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 2] += pâ * face.nÌ[2] * face.ÎS / cells[face.owner].Î©
        âÎpâx0[face.owner, 3] += pâ * face.nÌ[3] * face.ÎS / cells[face.owner].Î©
    end

    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal

        Ïâ = cells[face.owner].var[ð.Ï]
        Ïáµ£ = cells[face.neighbour].var[ð.Ï]
        pâ = cells[face.owner].var[ð.p]
        páµ£ = cells[face.neighbour].var[ð.p]
        uâ = cells[face.owner].var[ð.u]
        uáµ£ = cells[face.neighbour].var[ð.u]
        vâ = cells[face.owner].var[ð.v]
        váµ£ = cells[face.neighbour].var[ð.v]
        Î¼â = cells[face.owner].var[ð.Î¼]
        Î¼áµ£ = cells[face.neighbour].var[ð.Î¼]
        âÏâTâ = cells[face.owner].var[ð.âÏâT]
        âÏâTáµ£ = cells[face.neighbour].var[ð.âÏâT]
        âHââTâ = cells[face.owner].var[ð.âHââT]
        âHââTáµ£ = cells[face.neighbour].var[ð.âHââT]
        Hââ = cells[face.owner].var[ð.Hâ]
        Hâáµ£ = cells[face.neighbour].var[ð.Hâ]
        Uââ = uâ * face.nÌ[1] + vâ * face.nÌ[2]
        Uâáµ£ = uáµ£ * face.nÌ[1] + váµ£ * face.nÌ[2]
        Uâ = 0.5 * (Uââ + Uâáµ£)
        ÎS = face.ÎS

        centerâ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centeráµ£ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ÎLR = norm(centeráµ£ - centerâ)

        invÏÎt = 0.5 * (1.0/Ïâ + 1.0/Ïáµ£) * ð.Ît
        
        # Rhie-Chow
        Uâ += 0.5 * ð.Ît / Ïâ * âÎpâx0[face.owner, 1] * face.nÌ[1]
        Uâ += 0.5 * ð.Ît / Ïâ * âÎpâx0[face.owner, 2] * face.nÌ[2]
        Uâ += 0.5 * ð.Ît / Ïâ * âÎpâx0[face.owner, 3] * face.nÌ[3]
        Uâ += 0.5 * ð.Ît / Ïáµ£ * âÎpâx0[face.neighbour, 1] * face.nÌ[1]
        Uâ += 0.5 * ð.Ît / Ïáµ£ * âÎpâx0[face.neighbour, 2] * face.nÌ[2]
        Uâ += 0.5 * ð.Ît / Ïáµ£ * âÎpâx0[face.neighbour, 3] * face.nÌ[3]
        Uâ -= invÏÎt * (páµ£-pâ) / ÎLR

        Wâ = 0.5 * (1.0 + sign(Uâ))
        Wáµ£ = 1.0 - Wâ

        Ïâ = Wâ * Ïâ + Wáµ£ * Ïáµ£
        uâ = Wâ * uâ + Wáµ£ * uáµ£
        vâ = Wâ * vâ + Wáµ£ * váµ£
        Î¼â = Wâ * Î¼â + Wáµ£ * Î¼áµ£
#=        
        tmp_A_var_convL = Wâ * âÏâTâ * Hââ * Uâ * ÎS
        tmp_A_var_convR = Wáµ£ * âÏâTáµ£ * Hâáµ£ * Uâ * ÎS
        
        tmp_A_var_convL += Wâ * Ïâ * âHââTâ * Uâ * ÎS
        tmp_A_var_convR += Wáµ£ * Ïáµ£ * âHââTáµ£ * Uâ * ÎS

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
        covflux = ( Wâ * Ïâ * Hââ + Wáµ£ * Ïáµ£ * Hâáµ£ ) * Uâ * ÎS
        B[face.owner] -= covflux
        B[face.neighbour] += covflux

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        Ïâ = cells[face.owner].var[ð.Ï]
        ÎS = face.ÎS

        Uâ = 0.0
        Uâ += cells[face.owner].var[ð.u]*face.nÌ[1]
        Uâ += cells[face.owner].var[ð.v]*face.nÌ[2]
        Uâ += cells[face.owner].var[ð.w]*face.nÌ[3]
        Uâ0 = Uâ

        invU = cells[face.owner].var[ð.u] - Uâ * face.nÌ[1]
        invV = cells[face.owner].var[ð.v] - Uâ * face.nÌ[2]
        invW = cells[face.owner].var[ð.w] - Uâ * face.nÌ[3]

        Uâ = invU * face.nÌ[1]
        Uâ += invV * face.nÌ[2]
        Uâ += invW * face.nÌ[3]
        
        #A_vals[face.owner] += Ïâ * Uâ * ÎS

        Uâ = 0.0

        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    ps = MKLPardisoSolver()
    ÎT = solve(ps, A, B)
    




    relax = 0.5




    diagon = 1
    maximum_T = -1.e12
    for cell in cells

        cell.var[ð.T] += relax*ÎT[diagon]

        println(cell.var[ð.T])

        maximum_T = max(maximum_T,cell.var[ð.T])

        diagon += 1
    end


    return log10(norm(ÎT)/(maximum_T+1.e-20))
   

end
