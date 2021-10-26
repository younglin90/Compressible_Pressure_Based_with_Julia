function construct_A_matrix_explicit!(
    👉::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell}
)

    for i in 1:length(cells)

        ρ = cells[i].var[👉.ρ]
        u = cells[i].var[👉.u]
        v = cells[i].var[👉.v]
        w = cells[i].var[👉.w]
        T = cells[i].var[👉.T]
        Y₁ = cells[i].var[👉.Y₁]
        Hₜ = cells[i].var[👉.Hₜ]
        ∂ρ∂p = cells[i].var[👉.∂ρ∂p]
        ∂ρ∂T = cells[i].var[👉.∂ρ∂T]
        ∂Hₜ∂p = cells[i].var[👉.∂Hₜ∂p]
        ∂Hₜ∂T = cells[i].var[👉.∂Hₜ∂T]
        ∂ρ∂Y₁ = cells[i].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ = cells[i].var[👉.∂Hₜ∂Y₁]

        T = zeros(Float64, 6, 6)
        T[1, 1] = ∂ρ∂p
        T[2, 1] = ∂ρ∂p*u
        T[3, 1] = ∂ρ∂p*v
        T[4, 1] = ∂ρ∂p*w
        T[5, 1] = ∂ρ∂p*Hₜ + ρ*∂Hₜ∂p - 1.0
        
        T[2, 2] = ρ
        T[5, 2] = ρ*u
        
        T[3, 3] = ρ
        T[5, 3] = ρ*v
        
        T[4, 4] = ρ
        T[5, 4] = ρ*w
        
        T[1, 5] = ∂ρ∂T
        T[2, 5] = ∂ρ∂T*u
        T[3, 5] = ∂ρ∂T*v
        T[4, 5] = ∂ρ∂T*w
        T[5, 5] = ∂ρ∂T*Hₜ + ρ*∂Hₜ∂T

        T[6, 1] = ∂ρ∂p*Y₁
        T[6, 5] = ∂ρ∂T*Y₁

        T[6, 6] = ∂ρ∂Y₁*Y₁
        T[6, 6] += ρ

        T[1, 6] = ∂ρ∂Y₁
        T[2, 6] = ∂ρ∂Y₁*u
        T[3, 6] = ∂ρ∂Y₁*v
        T[4, 6] = ∂ρ∂Y₁*w
        T[5, 6] = ∂ρ∂Y₁*Hₜ + ρ*∂Hₜ∂Y₁

        P = deepcopy(T)
        β = 1.0/cells[i].var[👉.Vᵣ]^2.0
            + ∂ρ∂p - 1.0/cells[i].var[👉.c]^2.0
        P[1, 1] = β
        P[2, 1] = β*u
        P[3, 1] = β*v
        P[4, 1] = β*w
        P[5, 1] = β*Hₜ + ρ*∂Hₜ∂p - 1.0
        P[6, 1] = β*Y₁

        A[i, :, :] = T[:, :]./👉.Δt


    end


end



function linear_solver_explicit!(
    A::Array{Float64},
    x::Array{Float64},
    B::Array{Float64}
)




#@time begin
    n = length(A[:, 1, 1])
    for i in 1:n
        #invA = inv(A[i, :, :])
        #x[i, :] = invA[:, :]*B[i, :]
        
        x[i, :] = A[i, :, :]\B[i, :]

    end
#end

end



function flux!(
    👉::controls, 
    RHS::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
)


    for face in faces_internal
        cpi2 = 0.0
		w₁ = 0.0
		pₗ = face.varₗ[👉.p]
		ρₗ = face.varₗ[👉.ρ]
		cₗ = face.varₗ[👉.c]
		pᵣ = face.varᵣ[👉.p]
		ρᵣ = face.varᵣ[👉.ρ]
		cᵣ = face.varᵣ[👉.c]
		preLs = abs(pₗ) + 0.1 * ρₗ*cₗ*cₗ
		preRs = abs(pᵣ) + 0.1 * ρᵣ*cᵣ*cᵣ
		w₂ = min(preLs/preRs,preRs/preLs)

        Uₙ = 0.0

        flux = zeros(Float64, 6, 1)
        flux = KT_KNP(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            w₁, w₂, cpi2, Uₙ,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )


        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω #* 👉.Δt
        RHS[face.neighbour, :] += flux[:]*face.ΔS/cells[face.neighbour].Ω #* 👉.Δt

    end


    for face in faces_boundary

		Uₙ = 0.0

        flux = zeros(Float64, 6, 1)
        flux = KT_KNP(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            1.0, 1.0, 1.0, Uₙ,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )

        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω #* 👉.Δt

    end

	

end






function KT_KNP(
    pₗ,uₗ,vₗ,wₗ,Tₗ,Y₁ₗ,ρₗ,Hₜₗ,cₗ,
    pᵣ,uᵣ,vᵣ,wᵣ,Tᵣ,Y₁ᵣ,ρᵣ,Hₜᵣ,cᵣ,
    Lco,Uco,Δt,
    w₁,w₂,cpi2, Uₙ,
    nx,ny,nz
)

	Uₙₗ = uₗ*nx + vₗ*ny + wₗ*nz
	Uₙᵣ = uᵣ*nx + vᵣ*ny + wᵣ*nz

	Fmax =  max(max(Uₙₗ+cₗ,Uₙᵣ+cᵣ),0.0)
	Fmin = -min(min(Uₙₗ-cₗ,Uₙᵣ-cᵣ),0.0)

	α⁺ = Fmax / (Fmax+Fmin)
	α⁻ = Fmin / (Fmax+Fmin)
	α⁺⁻ = Fmax*Fmin / (Fmax+Fmin)

	ϕ⁺ = ρₗ * (α⁺ * Uₙₗ + α⁺⁻)
	ϕ⁻ = ρᵣ * (α⁻ * Uₙᵣ - α⁺⁻)

	ṁₗ = ϕ⁺
	ṁᵣ = ϕ⁻
	
	#ṁ = ϕ⁺ + ϕ⁻
	#ṁₗ = 0.5*(ṁ+abs(ṁ))
	#ṁᵣ = 0.5*(ṁ-abs(ṁ))

	
	pₗᵣ = 0.5*(pₗ + pᵣ)

	# comp. convective flux
	flux = zeros(Float64,6,1)
	flux[1] = ṁₗ + ṁᵣ
	flux[2] = ṁₗ*uₗ + ṁᵣ*uᵣ + pₗᵣ*nx
	flux[3] = ṁₗ*vₗ + ṁᵣ*vᵣ + pₗᵣ*ny
	flux[4] = ṁₗ*wₗ + ṁᵣ*wᵣ + pₗᵣ*nz
	flux[5] = ṁₗ*Hₜₗ + ṁᵣ*Hₜᵣ
	flux[6] = ṁₗ*Y₁ₗ + ṁᵣ*Y₁ᵣ

	return flux


end








