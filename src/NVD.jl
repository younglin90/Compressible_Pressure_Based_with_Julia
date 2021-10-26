
function calcUpwind!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            γ = 0.0
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end
	

function calcCentral!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            γ = 0.5
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end


function calcQUICK!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            C̃ = (αD-αU)/(abs(αA-αU)+1.e-200)*(αA>αU ? 1.0 : -1.0)
            γ = 0.0
            if 0.0 <= C̃ < 1.0
                C̃f = 0.375 + 0.75 * C̃
                γ = (C̃f-C̃)/(1.0-C̃)
            end
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end
	

function calcMinmod!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            C̃ = (αD-αU)/(abs(αA-αU)+1.e-200)*(αA>αU ? 1.0 : -1.0)
            γ = 0.0
            if 0.0 <= C̃ < 0.5
                C̃f = 1.5 * C̃
                γ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.5 <= C̃ < 1.0
                C̃f = 0.5 * C̃ + 0.5
                γ = (C̃f-C̃)/(1.0-C̃)
            end
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end
	



function calcVanLeer!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            C̃ = (αD-αU)/(abs(αA-αU)+1.e-200)*(αA>αU ? 1.0 : -1.0)
            γ = 0.0
            if 0.0 <= C̃ < 1.0
                C̃f = C̃ + (C̃*(1.0-C̃))/(2.0-4.0*C̃+4.0*C̃^2)
                γ = (C̃f-C̃)/(1.0-C̃)
            end
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end
	





function calcSUPERBEE!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
            end

            C̃ = (αD-αU)/(abs(αA-αU)+1.e-200)*(αA>αU ? 1.0 : -1.0)
            γ = 0.0
            if 0.0 <= C̃ < 0.5
                C̃f = 0.5 * C̃ + 0.5
                γ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.5 <= C̃ < 0.6666
                C̃f = 1.5 * C̃
                γ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.6666 <= C̃ < 1.0
                C̃f = 1.0
                γ = (C̃f-C̃)/(1.0-C̃)
            end
            
            if iter==1
                face.varₗ[ϕ] = γ*αA + (1.0-γ)*αD
            else
                face.varᵣ[ϕ] = γ*αA + (1.0-γ)*αD
            end
        end

    end


end
	







function calcMSTACS!(
    👉::controls,
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    corant = zeros(Float64, length(cells))
    for face in faces_internal
        corantNo_L = 👉.Δt * face.ΔS / cells[face.owner].Ω * 
            sqrt(cells[face.owner].var[👉.u]^2 + cells[face.owner].var[👉.v]^2 + cells[face.owner].var[👉.w]^2 + 1.e-100)
        corantNo_R = 👉.Δt * face.ΔS / cells[face.neighbour].Ω * 
            sqrt(cells[face.neighbour].var[👉.u]^2 + cells[face.neighbour].var[👉.v]^2 + cells[face.neighbour].var[👉.w]^2 + 1.e-100)
        corantNo_max = max(corantNo_L,corantNo_R)
        corant[face.owner] = max(corant[face.owner],corantNo_max)
        corant[face.neighbour] = max(corant[face.neighbour],corantNo_max)
    end



    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ₋₁ = ϕᵢ
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₋₁ -= ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₂ = ϕᵢ₊₁
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₂ += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)

        for iter in 1:2
            tmp1 = 0.0
            cosθ = 0.0
            dist = (cells[face.neighbour].x - cells[face.owner].x)^2
            dist += (cells[face.neighbour].y - cells[face.owner].y)^2
            dist += (cells[face.neighbour].z - cells[face.owner].z)^2
            if iter==1
                αU = ϕᵢ₋₁; αD = ϕᵢ; αA = ϕᵢ₊₁
                corantNo = corant[face.owner]
                cosθ = ∂ϕ∂x[face.owner, 1]^2 + ∂ϕ∂x[face.owner, 2]^2 + ∂ϕ∂x[face.owner, 3]^2
                cosθ = sqrt(cosθ)*sqrt(dist)
    
                tmp1 = ∂ϕ∂x[face.owner, 1] * (cells[face.neighbour].x - cells[face.owner].x)
                tmp1 += ∂ϕ∂x[face.owner, 2] * (cells[face.neighbour].y - cells[face.owner].y)
                tmp1 += ∂ϕ∂x[face.owner, 3] * (cells[face.neighbour].z - cells[face.owner].z)
            else
                αU = ϕᵢ₊₂; αD = ϕᵢ₊₁; αA = ϕᵢ
                corantNo = corant[face.owner]
                cosθ = ∂ϕ∂x[face.neighbour, 1]^2 + ∂ϕ∂x[face.neighbour, 2]^2 + ∂ϕ∂x[face.neighbour, 3]^2
                cosθ = sqrt(cosθ)*sqrt(dist)
    
                tmp1 = ∂ϕ∂x[face.neighbour, 1] * (cells[face.neighbour].x - cells[face.owner].x)
                tmp1 += ∂ϕ∂x[face.neighbour, 2] * (cells[face.neighbour].y - cells[face.owner].y)
                tmp1 += ∂ϕ∂x[face.neighbour, 3] * (cells[face.neighbour].z - cells[face.owner].z)
            end

            if cosθ > 0.0
                tmp1 = abs(tmp1) / cosθ
            end

            γ = min(cosθ^4,1.0)
            
            C̃ = (αD-αU)/(abs(αA-αU)+1.e-200)*(αA>αU ? 1.0 : -1.0)
            γ₁ = 0.0
            if 0.0 <= C̃ < 1.0 && 0.0 <= corantNo <= 0.33
                C̃f = min(1.0, C̃/corantNo)
                γ₁ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.0 <= C̃ < 1.0 && 0.33 < corantNo <= 1.0
                C̃f = min(1.0,3.0*C̃)
                γ₁ = (C̃f-C̃)/(1.0-C̃)
            else
                γ₁ = 0.0
            end
            
            vfCompressive = γ₁*αA + (1.0-γ₁)*αD

            γ₂ = 0.0
            if 0.0 <= C̃ < 0.2
                C̃f = min(1.0, 3.0*C̃)
                γ₂ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.2 <= C̃ < 0.5
                C̃f = 0.5*C̃ + 0.5
                γ₂ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.5 <= C̃ < 0.8333
                C̃f = 0.75*C̃ + 0.375
                γ₂ = (C̃f-C̃)/(1.0-C̃)
            elseif 0.8333 <= C̃ < 1.0
                C̃f = 1.0
                γ₂ = (C̃f-C̃)/(1.0-C̃)
            else
                γ₂ = 0.0
            end

            vfDiffusive = γ₂*αA + (1.0-γ₂)*αD

            if iter==1
                face.varₗ[ϕ] = γ*vfCompressive + (1.0-γ)*vfDiffusive
            else
                face.varᵣ[ϕ] = γ*vfCompressive + (1.0-γ)*vfDiffusive
            end
        end

    end


end
	









function calcLinearInterpolation!(
    ϕ::UInt32,
    ∂ϕ∂x::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face}
    )

    
    for face in faces_internal

        ϕᵢ = cells[face.owner].var[ϕ]
        ϕᵢ₊₁ = cells[face.neighbour].var[ϕ]

        ϕᵢ += ∂ϕ∂x[face.owner, 1] * (face.x - cells[face.owner].x)
        ϕᵢ += ∂ϕ∂x[face.owner, 2] * (face.y - cells[face.owner].y)
        ϕᵢ += ∂ϕ∂x[face.owner, 3] * (face.z - cells[face.owner].z)

        ϕᵢ₊₁ += ∂ϕ∂x[face.neighbour, 1] * (face.x - cells[face.neighbour].x)
        ϕᵢ₊₁ += ∂ϕ∂x[face.neighbour, 2] * (face.y - cells[face.neighbour].y)
        ϕᵢ₊₁ += ∂ϕ∂x[face.neighbour, 3] * (face.z - cells[face.neighbour].z)
        
        #=
        ϕᵢ += ∂ϕ∂x[face.owner, 1] * 0.5 * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ += ∂ϕ∂x[face.owner, 2] * 0.5 * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ += ∂ϕ∂x[face.owner, 3] * 0.5 * (cells[face.neighbour].z - cells[face.owner].z)

        ϕᵢ₊₁ -= ∂ϕ∂x[face.neighbour, 1] * 0.5 * (cells[face.neighbour].x - cells[face.owner].x)
        ϕᵢ₊₁ -= ∂ϕ∂x[face.neighbour, 2] * 0.5 * (cells[face.neighbour].y - cells[face.owner].y)
        ϕᵢ₊₁ -= ∂ϕ∂x[face.neighbour, 3] * 0.5 * (cells[face.neighbour].z - cells[face.owner].z)
        =#
        
        face.varₗ[ϕ] = ϕᵢ
        face.varᵣ[ϕ] = ϕᵢ₊₁

    end


end
	


