
function reconstruction!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
    )



    ∂p∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.p] + cells[face.neighbour].var[👉.p])
        ∂p∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂p∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂p∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂p∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂p∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂p∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.p]
        ∂p∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂p∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂p∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_p == "upwind"
        calcUpwind!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "central"
        calcCentral!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "quick"
        calcQUICK!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "minmod"
        calcMinmod!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "vanleer"
        calcVanLeer!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "superbee"
        calcSUPERBEE!(👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "mstacs"
        calcMSTACS!(👉, 👉.p, ∂p∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_p == "linear"
        calcLinearInterpolation!(👉.p, ∂p∂x, cells, faces_internal)
    end
    
    ∂u∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.u] + cells[face.neighbour].var[👉.u])
        ∂u∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂u∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂u∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂u∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂u∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂u∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.u]
        ∂u∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂u∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂u∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_U == "upwind"
        calcUpwind!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "central"
        calcCentral!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "quick"
        calcQUICK!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "minmod"
        calcMinmod!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "vanleer"
        calcVanLeer!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "superbee"
        calcSUPERBEE!(👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "mstacs"
        calcMSTACS!(👉, 👉.u, ∂u∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "linear"
        calcLinearInterpolation!(👉.u, ∂u∂x, cells, faces_internal)
    end
    
    ∂v∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.v] + cells[face.neighbour].var[👉.v])
        ∂v∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂v∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂v∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂v∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂v∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂v∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.v]
        ∂v∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂v∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂v∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_U == "upwind"
        calcUpwind!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "central"
        calcCentral!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "quick"
        calcQUICK!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "minmod"
        calcMinmod!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "vanleer"
        calcVanLeer!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "superbee"
        calcSUPERBEE!(👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "mstacs"
        calcMSTACS!(👉, 👉.v, ∂v∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_U == "linear"
        calcLinearInterpolation!(👉.v, ∂v∂x, cells, faces_internal)
    end

    

    ∂T∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.T] + cells[face.neighbour].var[👉.T])
        ∂T∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂T∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂T∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂T∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂T∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂T∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.T]
        ∂T∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂T∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂T∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_T == "upwind"
        calcUpwind!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "central"
        calcCentral!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "quick"
        calcQUICK!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "minmod"
        calcMinmod!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "vanleer"
        calcVanLeer!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "superbee"
        calcSUPERBEE!(👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "mstacs"
        calcMSTACS!(👉, 👉.T, ∂T∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_T == "linear"
        calcLinearInterpolation!(👉.T, ∂T∂x, cells, faces_internal)
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
    elseif 👉.spatial_discretizationScheme_Y == "mstacs"
        calcMSTACS!(👉, 👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "linear"
        calcLinearInterpolation!(👉.Y₁, ∂Y₁∂x, cells, faces_internal)
    end


    for face in faces_internal
        pₙ = face.varₗ[👉.p]
        uₙ = face.varₗ[👉.u]
        vₙ = face.varₗ[👉.v]
        wₙ = 0.0
        Tₙ = face.varₗ[👉.T]
        Y₁ₙ = face.varₗ[👉.Y₁]
        ρₙ, Hₜₙ, cₙ, ∂ρ∂pₙ, ∂ρ∂Tₙ, ∂ρ∂Y₁ₙ, ∂Hₜ∂Y₁ₙ = faceEOS!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)
        face.varₗ[👉.ρ] = ρₙ
        face.varₗ[👉.Hₜ] = Hₜₙ
        face.varₗ[👉.c] = cₙ
        face.varₗ[👉.∂ρ∂p] = ∂ρ∂pₙ
        face.varₗ[👉.∂ρ∂T] = ∂ρ∂Tₙ
        face.varₗ[👉.∂ρ∂Y₁] = ∂ρ∂Y₁ₙ
        face.varₗ[👉.∂Hₜ∂Y₁] = ∂Hₜ∂Y₁ₙ
        
        pₙ = face.varᵣ[👉.p]
        uₙ = face.varᵣ[👉.u]
        vₙ = face.varᵣ[👉.v]
        wₙ = 0.0
        Tₙ = face.varᵣ[👉.T]
        Y₁ₙ = face.varᵣ[👉.Y₁]
        ρₙ, Hₜₙ, cₙ, ∂ρ∂pₙ, ∂ρ∂Tₙ, ∂ρ∂Y₁ₙ, ∂Hₜ∂Y₁ₙ = faceEOS!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)
        face.varᵣ[👉.ρ] = ρₙ
        face.varᵣ[👉.Hₜ] = Hₜₙ
        face.varᵣ[👉.c] = cₙ
        face.varᵣ[👉.∂ρ∂p] = ∂ρ∂pₙ
        face.varᵣ[👉.∂ρ∂T] = ∂ρ∂Tₙ
        face.varᵣ[👉.∂ρ∂Y₁] = ∂ρ∂Y₁ₙ
        face.varᵣ[👉.∂Hₜ∂Y₁] = ∂Hₜ∂Y₁ₙ
    end
    #=

    ∂α₁∂x = zeros(Float64, length(cells), 3)
    for face in faces_internal
        ϕₙ = 0.5 * (cells[face.owner].var[👉.α₁] + cells[face.neighbour].var[👉.α₁])
        ∂α₁∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂α₁∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂α₁∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂α₁∂x[face.neighbour, 1] -= ϕₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂α₁∂x[face.neighbour, 2] -= ϕₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂α₁∂x[face.neighbour, 3] -= ϕₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end
    for face in faces_boundary
        ϕₙ = cells[face.owner].var[👉.α₁]
        ∂α₁∂x[face.owner, 1] += ϕₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂α₁∂x[face.owner, 2] += ϕₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂α₁∂x[face.owner, 3] += ϕₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end
    if 👉.spatial_discretizationScheme_Y == "upwind"
        calcUpwind!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "central"
        calcCentral!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "quick"
        calcQUICK!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "minmod"
        calcMinmod!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "vanleer"
        calcVanLeer!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "superbee"
        calcSUPERBEE!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "mstacs"
        calcMSTACS!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    elseif 👉.spatial_discretizationScheme_Y == "linear"
        calcLinearInterpolation!(👉.α₁, ∂α₁∂x, cells, faces_internal)
    end


    for face in faces_internal
        pₙ = face.varₗ[👉.p]
        uₙ = face.varₗ[👉.u]
        vₙ = face.varₗ[👉.v]
        wₙ = 0.0
        Tₙ = face.varₗ[👉.T]
        α₁ₙ = face.varₗ[👉.α₁]
        ρₙ, Hₜₙ, cₙ, ∂ρ∂pₙ, ∂ρ∂Tₙ, ∂ρ∂Y₁ₙ, ∂Hₜ∂Y₁ₙ = faceEOS_vf!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)
        face.varₗ[👉.ρ] = ρₙ
        face.varₗ[👉.Hₜ] = Hₜₙ
        face.varₗ[👉.c] = cₙ
        face.varₗ[👉.∂ρ∂p] = ∂ρ∂pₙ
        face.varₗ[👉.∂ρ∂T] = ∂ρ∂Tₙ
        face.varₗ[👉.∂ρ∂Y₁] = ∂ρ∂Y₁ₙ
        face.varₗ[👉.∂Hₜ∂Y₁] = ∂Hₜ∂Y₁ₙ
        
        pₙ = face.varᵣ[👉.p]
        uₙ = face.varᵣ[👉.u]
        vₙ = face.varᵣ[👉.v]
        wₙ = 0.0
        Tₙ = face.varᵣ[👉.T]
        α₁ₙ = face.varᵣ[👉.α₁]
        ρₙ, Hₜₙ, cₙ, ∂ρ∂pₙ, ∂ρ∂Tₙ, ∂ρ∂Y₁ₙ, ∂Hₜ∂Y₁ₙ = faceEOS_vf!(👉,pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)
        face.varᵣ[👉.ρ] = ρₙ
        face.varᵣ[👉.Hₜ] = Hₜₙ
        face.varᵣ[👉.c] = cₙ
        face.varᵣ[👉.∂ρ∂p] = ∂ρ∂pₙ
        face.varᵣ[👉.∂ρ∂T] = ∂ρ∂Tₙ
        face.varᵣ[👉.∂ρ∂Y₁] = ∂ρ∂Y₁ₙ
        face.varᵣ[👉.∂Hₜ∂Y₁] = ∂Hₜ∂Y₁ₙ
    end
=#



end