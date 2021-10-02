#= include("./module.jl")
include("./controls.jl")
using .mesh =#

function EOS!(
    👉::controls, cell::Vector{mesh.Cell}
)


    for i in cell
        ρᵢ = zeros(Float64, 2, 1)
        ρᵢ = zeros(Float64, 2, 1)
        Hₜᵢ = zeros(Float64, 2, 1)
        cᵢ = zeros(Float64, 2, 1)
        ∂ρ∂pᵢ = zeros(Float64, 2, 1)
        ∂ρ∂Tᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)

        i.var[👉.Y₁] = max(min(i.var[👉.Y₁],1.0),0.0)
        i.var[👉.Y₂] = 1.0 - i.var[👉.Y₁]

        j = 1
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosNASG(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])
        j = 2
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosIdeal(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])

        ρ = 1.0/(i.var[👉.Y₁]/ρᵢ[1]+i.var[👉.Y₂]/ρᵢ[2])

        i.var[👉.α₁] = ρ*i.var[👉.Y₁]/ρᵢ[1]
        i.var[👉.α₁] = max(min(i.var[👉.α₁],1.0),0.0)
        i.var[👉.α₂] = 1.0 - i.var[👉.α₁]

        i.var[👉.ρ] = ρ
        i.var[👉.Hₜ] = i.var[👉.Y₁]*Hₜᵢ[1] + i.var[👉.Y₂]*Hₜᵢ[2]

        #i.var[👉.∂ρ∂p] = i.var[👉.α₁]*∂ρ∂pᵢ[1] + i.var[👉.α₂]*∂ρ∂pᵢ[2]
        #i.var[👉.∂ρ∂T] = i.var[👉.α₁]*∂ρ∂Tᵢ[1] + i.var[👉.α₂]*∂ρ∂Tᵢ[2]
        i.var[👉.∂ρ∂p] = 
			ρ^2*(i.var[👉.Y₁]/ρᵢ[1]^2*∂ρ∂pᵢ[1] + 
				 i.var[👉.Y₂]/ρᵢ[2]^2*∂ρ∂pᵢ[2])
        i.var[👉.∂ρ∂T] = 
			ρ^2*(i.var[👉.Y₁]/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + 
				 i.var[👉.Y₂]/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

        i.var[👉.∂Hₜ∂p] = i.var[👉.Y₁]*∂Hₜ∂pᵢ[1] + i.var[👉.Y₂]*∂Hₜ∂pᵢ[2]
        i.var[👉.∂Hₜ∂T] = i.var[👉.Y₁]*∂Hₜ∂Tᵢ[1] + i.var[👉.Y₂]*∂Hₜ∂Tᵢ[2]

        i.var[👉.∂ρ∂Y₁] = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
        i.var[👉.∂Hₜ∂Y₁] = Hₜᵢ[1]-Hₜᵢ[2]

        i.var[👉.c] = i.var[👉.∂ρ∂p] + 1.0/ρ*i.var[👉.∂ρ∂T]/i.var[👉.∂Hₜ∂T]*(1.0-ρ*i.var[👉.∂Hₜ∂p])
        i.var[👉.c] = √(1.0/i.var[👉.c])


    end


end


function EOS_vf!(
    👉::controls, cell::Vector{mesh.Cell}
)


    for i in cell
        ρᵢ = zeros(Float64, 2, 1)
        ρᵢ = zeros(Float64, 2, 1)
        Hₜᵢ = zeros(Float64, 2, 1)
        cᵢ = zeros(Float64, 2, 1)
        ∂ρ∂pᵢ = zeros(Float64, 2, 1)
        ∂ρ∂Tᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)

        i.var[👉.α₁] = max(min(i.var[👉.α₁],1.0),0.0)
        i.var[👉.α₂] = 1.0 - i.var[👉.α₁]

		α₁ = i.var[👉.α₁]
		α₂ = i.var[👉.α₂]

        j = 1
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosNASG(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])
        j = 2
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosIdeal(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])

        ρ = α₁ * ρᵢ[1] + α₂ * ρᵢ[2]

		Y₁ = α₁*ρᵢ[1]/ρ
		Y₂= α₂*ρᵢ[2]/ρ

        i.var[👉.Y₁] = ρᵢ[1]*i.var[👉.α₁]/ρ
        i.var[👉.Y₁] = max(min(i.var[👉.Y₁],1.0),0.0)
        i.var[👉.Y₂] = 1.0 - i.var[👉.Y₁]

        i.var[👉.ρ] = ρ
		
        i.var[👉.Hₜ] = (α₁ * ρᵢ[1] * Hₜᵢ[1] + α₂ * ρᵢ[2] * Hₜᵢ[2]) / ρ

		∂ρ∂p = α₁ *∂ρ∂pᵢ[1] + α₂ * ∂ρ∂pᵢ[2]
        i.var[👉.∂ρ∂p] = ∂ρ∂p

		∂ρ∂T = α₁ *∂ρ∂Tᵢ[1] + α₂ * ∂ρ∂Tᵢ[2]
        i.var[👉.∂ρ∂T] = ∂ρ∂T

		i.var[👉.∂Hₜ∂p] = ( 
			1.0/ρ*(∂ρ∂pᵢ[1]*α₁-∂ρ∂p*Y₁)*Hₜᵢ[1] +
			1.0/ρ*(∂ρ∂pᵢ[2]*α₂-∂ρ∂p*Y₂)*Hₜᵢ[2] +
			Y₁*∂Hₜ∂pᵢ[1] + Y₂*∂Hₜ∂pᵢ[2]
		)
		#i.var[👉.∂Hₜ∂p] += ρᵢ[1] * α₁ * ∂Hₜ∂pᵢ[1] + ρᵢ[2] * α₂ * ∂Hₜ∂pᵢ[2]
		#i.var[👉.∂Hₜ∂p] -= ∂ρ∂p
		#i.var[👉.∂Hₜ∂p] /= ρ

		i.var[👉.∂Hₜ∂T] = ( 
			1.0/ρ*(∂ρ∂Tᵢ[1]*α₁-∂ρ∂T*Y₁)*Hₜᵢ[1] +
			1.0/ρ*(∂ρ∂Tᵢ[2]*α₂-∂ρ∂T*Y₂)*Hₜᵢ[2] +
			Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]
		)

		#i.var[👉.∂Hₜ∂T] = ∂ρ∂Tᵢ[1] * α₁ * Hₜᵢ[1] + ∂ρ∂Tᵢ[2] * α₂ * Hₜᵢ[2]
		#i.var[👉.∂Hₜ∂T] += ρᵢ[1] * α₁ * ∂Hₜ∂Tᵢ[1] + ρᵢ[2] * α₂ * ∂Hₜ∂Tᵢ[2]
		#i.var[👉.∂Hₜ∂T] -= ∂ρ∂T
		#i.var[👉.∂Hₜ∂T] /= ρ

        #i.var[👉.∂ρ∂Y₁] = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
        #i.var[👉.∂Hₜ∂Y₁] = Hₜᵢ[1]-Hₜᵢ[2]
        #i.var[👉.∂ρ∂Y₁] = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
        #i.var[👉.∂Hₜ∂Y₁] = Hₜᵢ[1]-Hₜᵢ[2]
        i.var[👉.∂ρ∂Y₁] = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
        i.var[👉.∂Hₜ∂Y₁] = Hₜᵢ[1]-Hₜᵢ[2]

        i.var[👉.c] = i.var[👉.∂ρ∂p] + 1.0/ρ*i.var[👉.∂ρ∂T]/i.var[👉.∂Hₜ∂T]*(1.0-ρ*i.var[👉.∂Hₜ∂p])
        i.var[👉.c] = √(1.0/i.var[👉.c])


    end


end



function faceEOS!(
	p::Float64,u::Float64,v::Float64,w::Float64,T::Float64,Y₁::Float64
)


	ρᵢ = zeros(Float64, 2, 1)
	ρᵢ = zeros(Float64, 2, 1)
	Hₜᵢ = zeros(Float64, 2, 1)
	cᵢ = zeros(Float64, 2, 1)
	∂ρ∂pᵢ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)

	Y₁ = max(min(Y₁,1.0),0.0)
	Y₂ = 1.0 - Y₁
	#Y₁ = max(min(Y₁,1.0),0.0)
	#Y₂ = 1.0 - Y₁

	j = 1
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosNASG(p, u, v, w, T)
	j = 2
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosIdeal(p, u, v, w, T)

	ρ = 1.0/(Y₁/ρᵢ[1]+Y₂/ρᵢ[2])

	α₁ = ρ*Y₁/ρᵢ[1]
	α₁ = max(min(α₁,1.0),0.0)
	α₂ = 1.0 - α₁

	Hₜ = Y₁*Hₜᵢ[1] + Y₂*Hₜᵢ[2]

	#∂ρ∂p = α₁*∂ρ∂pᵢ[1] + α₂*∂ρ∂pᵢ[2]
	#∂ρ∂T = α₁*∂ρ∂Tᵢ[1] + α₂*∂ρ∂Tᵢ[2]
    ∂ρ∂p = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂pᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂pᵢ[2])
    ∂ρ∂T = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

	∂Hₜ∂p = Y₁*∂Hₜ∂pᵢ[1] + Y₂*∂Hₜ∂pᵢ[2]
	∂Hₜ∂T = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]

	#∂ρ∂Y₁ = -ρ*ρ*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
	#∂Hₜ∂Y₁ = Hₜᵢ[1]-Hₜᵢ[2]

	c = ∂ρ∂p + 1.0/ρ*∂ρ∂T/∂Hₜ∂T*(1.0-ρ*∂Hₜ∂p)
	c = √(1.0/c)

	return ρ, Hₜ, c


end



function faceEOS_vf!(
	p::Float64,u::Float64,v::Float64,w::Float64,T::Float64,α₁::Float64
)

	ρᵢ = zeros(Float64, 2, 1)
	ρᵢ = zeros(Float64, 2, 1)
	Hₜᵢ = zeros(Float64, 2, 1)
	cᵢ = zeros(Float64, 2, 1)
	∂ρ∂pᵢ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)




	α₁ = max(min(α₁,1.0),0.0)
	α₂ = 1.0 - α₁

	j = 1
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosNASG(p, u, v, w, T)
	j = 2
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosIdeal(p, u, v, w, T)


	

	ρ = α₁ * ρᵢ[1] + α₂ * ρᵢ[2]

	Y₁ = α₁*ρᵢ[1]/ρ
	Y₂= α₂*ρᵢ[2]/ρ
	
	Hₜ = (α₁ * ρᵢ[1] * Hₜᵢ[1] + α₂ * ρᵢ[2] * Hₜᵢ[2]) / ρ

	∂ρ∂p = α₁ *∂ρ∂pᵢ[1] + α₂ * ∂ρ∂pᵢ[2]

	∂ρ∂T = α₁ *∂ρ∂Tᵢ[1] + α₂ * ∂ρ∂Tᵢ[2]

	∂Hₜ∂p = ( 
		1.0/ρ*(∂ρ∂pᵢ[1]*α₁-∂ρ∂p*Y₁)*Hₜᵢ[1] +
		1.0/ρ*(∂ρ∂pᵢ[2]*α₂-∂ρ∂p*Y₂)*Hₜᵢ[2] +
		Y₁*∂Hₜ∂pᵢ[1] + Y₂*∂Hₜ∂pᵢ[2]
	)
	∂Hₜ∂T = ( 
		1.0/ρ*(∂ρ∂Tᵢ[1]*α₁-∂ρ∂T*Y₁)*Hₜᵢ[1] +
		1.0/ρ*(∂ρ∂Tᵢ[2]*α₂-∂ρ∂T*Y₂)*Hₜᵢ[2] +
		Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]
	)

	c = ∂ρ∂p + 1.0/ρ*∂ρ∂T/∂Hₜ∂T*(1.0-ρ*∂Hₜ∂p)
	c = √(1.0/c)

	return ρ, Hₜ, c


end




function EOS_ACID(
	pₗ::Float64,pᵣ::Float64,
	uₗ::Float64,uᵣ::Float64,
	vₗ::Float64,vᵣ::Float64,
	wₗ::Float64,wᵣ::Float64,
	Tₗ::Float64,Tᵣ::Float64,
	Y₁ₗ::Float64,Y₁ᵣ::Float64
)

	ρᵢₗ = zeros(Float64, 2, 1)
	ρᵢₗ = zeros(Float64, 2, 1)
	Hₜᵢₗ = zeros(Float64, 2, 1)
	cᵢₗ = zeros(Float64, 2, 1)
	∂ρ∂pᵢₗ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢₗ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢₗ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢₗ = zeros(Float64, 2, 1)
	
	ρᵢᵣ = zeros(Float64, 2, 1)
	ρᵢᵣ = zeros(Float64, 2, 1)
	Hₜᵢᵣ = zeros(Float64, 2, 1)
	cᵢᵣ = zeros(Float64, 2, 1)
	∂ρ∂pᵢᵣ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢᵣ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢᵣ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢᵣ = zeros(Float64, 2, 1)

	Y₂ₗ = 1.0 - Y₁ₗ
	Y₂ᵣ = 1.0 - Y₁ᵣ

	j = 1
	(ρᵢₗ[j], Hₜᵢₗ[j], cᵢₗ[j], ∂ρ∂pᵢₗ[j], ∂ρ∂Tᵢₗ[j], ∂Hₜ∂pᵢₗ[j], ∂Hₜ∂Tᵢₗ[j]) = 
	eosNASG(pₗ, uₗ, vₗ, wₗ, Tₗ)
	j = 2
	(ρᵢₗ[j], Hₜᵢₗ[j], cᵢₗ[j], ∂ρ∂pᵢₗ[j], ∂ρ∂Tᵢₗ[j], ∂Hₜ∂pᵢₗ[j], ∂Hₜ∂Tᵢₗ[j]) = 
	eosIdeal(pₗ, uₗ, vₗ, wₗ, Tₗ)
	
	j = 1
	(ρᵢᵣ[j], Hₜᵢᵣ[j], cᵢᵣ[j], ∂ρ∂pᵢᵣ[j], ∂ρ∂Tᵢᵣ[j], ∂Hₜ∂pᵢᵣ[j], ∂Hₜ∂Tᵢᵣ[j]) = 
	eosNASG(pᵣ, uᵣ, vᵣ, wᵣ, Tᵣ)
	j = 2
	(ρᵢᵣ[j], Hₜᵢᵣ[j], cᵢᵣ[j], ∂ρ∂pᵢᵣ[j], ∂ρ∂Tᵢᵣ[j], ∂Hₜ∂pᵢᵣ[j], ∂Hₜ∂Tᵢᵣ[j]) = 
	eosIdeal(pᵣ, uᵣ, vᵣ, wᵣ, Tᵣ)

	ρₗ_ACID = 1.0/(Y₁ᵣ/ρᵢₗ[1]+Y₂ᵣ/ρᵢₗ[2])
	ρᵣ_ACID = 1.0/(Y₁ₗ/ρᵢᵣ[1]+Y₂ₗ/ρᵢᵣ[2])
	
    ∂ρ∂pₗ_ACID = ρₗ_ACID^2*(Y₁ᵣ/ρᵢₗ[1]^2*∂ρ∂pᵢₗ[1] + Y₂ᵣ/ρᵢₗ[2]^2*∂ρ∂pᵢₗ[2])
    ∂ρ∂pᵣ_ACID = ρᵣ_ACID^2*(Y₁ₗ/ρᵢᵣ[1]^2*∂ρ∂pᵢᵣ[1] + Y₂ₗ/ρᵢᵣ[2]^2*∂ρ∂pᵢᵣ[2])

	∂Hₜ∂pₗ_ACID = Y₁ᵣ*∂Hₜ∂pᵢₗ[1] + Y₂ᵣ*∂Hₜ∂pᵢₗ[2]
	∂Hₜ∂pᵣ_ACID = Y₁ₗ*∂Hₜ∂pᵢᵣ[1] + Y₂ₗ*∂Hₜ∂pᵢᵣ[2]
	Hₜₗ_ACID = Y₁ᵣ*Hₜᵢₗ[1] + Y₂ᵣ*Hₜᵢₗ[2]
	Hₜᵣ_ACID = Y₁ₗ*Hₜᵢᵣ[1] + Y₂ₗ*Hₜᵢᵣ[2]
	∂ρ∂Tₗ_ACID = ρₗ_ACID^2*(Y₁ᵣ/ρᵢₗ[1]^2*∂ρ∂Tᵢₗ[1] + Y₂ᵣ/ρᵢₗ[2]^2*∂ρ∂Tᵢₗ[2])
	∂ρ∂Tᵣ_ACID = ρᵣ_ACID^2*(Y₁ₗ/ρᵢᵣ[1]^2*∂ρ∂Tᵢᵣ[1] + Y₂ₗ/ρᵢᵣ[2]^2*∂ρ∂Tᵢᵣ[2])
	∂Hₜ∂Tₗ_ACID = Y₁ᵣ*∂Hₜ∂Tᵢₗ[1] + Y₂ᵣ*∂Hₜ∂Tᵢₗ[2]
	∂Hₜ∂Tᵣ_ACID = Y₁ₗ*∂Hₜ∂Tᵢᵣ[1] + Y₂ₗ*∂Hₜ∂Tᵢᵣ[2]
	
	∂ρ∂Y₁ₗ_ACID = -ρₗ_ACID^2*(1.0/ρᵢₗ[1] + 1.0/ρᵢₗ[2])
	∂ρ∂Y₁ᵣ_ACID = -ρᵣ_ACID^2*(1.0/ρᵢᵣ[1] + 1.0/ρᵢᵣ[2])
	#∂Hₜ∂Y₁ₗ_ACID = Y₁ᵣ*∂Hₜ∂Tᵢᵣ[1] + Y₂ᵣ*∂Hₜ∂Tᵢᵣ[2]
	#∂Hₜ∂Y₁ᵣ_ACID = Y₁ₗ*∂Hₜ∂Tᵢᵣ[1] + Y₂ₗ*∂Hₜ∂Tᵢᵣ[2]

	#Hₜ = Y₁*Hₜᵢ[1] + Y₂*Hₜᵢ[2]

	#∂ρ∂p = α₁*∂ρ∂pᵢ[1] + α₂*∂ρ∂pᵢ[2]
	#∂ρ∂T = α₁*∂ρ∂Tᵢ[1] + α₂*∂ρ∂Tᵢ[2]
    #∂ρ∂pₗ_ACID = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂pᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂pᵢ[2])
    #∂ρ∂T = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

	#∂Hₜ∂p = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]
	#∂Hₜ∂T = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]

	#∂ρ∂Y₁ = -ρ*ρ*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
	#∂Hₜ∂Y₁ = Hₜᵢ[1]-Hₜᵢ[2]

	#c = ∂ρ∂p + 1.0/ρ*∂ρ∂T/∂Hₜ∂T*(1.0-ρ*∂Hₜ∂p)
	#c = √(1.0/c)

	return ρₗ_ACID,ρᵣ_ACID,∂ρ∂pₗ_ACID,∂ρ∂pᵣ_ACID, ∂Hₜ∂pₗ_ACID, ∂Hₜ∂pᵣ_ACID,
	Hₜₗ_ACID, Hₜᵣ_ACID, ∂ρ∂Tₗ_ACID, ∂ρ∂Tᵣ_ACID, ∂Hₜ∂Tₗ_ACID, ∂Hₜ∂Tᵣ_ACID,
	∂ρ∂Y₁ₗ_ACID, ∂ρ∂Y₁ᵣ_ACID


end





function EOS_ACID_vf(
	pₗ::Float64,pᵣ::Float64,
	uₗ::Float64,uᵣ::Float64,
	vₗ::Float64,vᵣ::Float64,
	wₗ::Float64,wᵣ::Float64,
	Tₗ::Float64,Tᵣ::Float64,
	α₁ₗ::Float64,α₁ᵣ::Float64
)

	ρᵢₗ = zeros(Float64, 2, 1)
	ρᵢₗ = zeros(Float64, 2, 1)
	Hₜᵢₗ = zeros(Float64, 2, 1)
	cᵢₗ = zeros(Float64, 2, 1)
	∂ρ∂pᵢₗ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢₗ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢₗ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢₗ = zeros(Float64, 2, 1)
	
	ρᵢᵣ = zeros(Float64, 2, 1)
	ρᵢᵣ = zeros(Float64, 2, 1)
	Hₜᵢᵣ = zeros(Float64, 2, 1)
	cᵢᵣ = zeros(Float64, 2, 1)
	∂ρ∂pᵢᵣ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢᵣ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢᵣ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢᵣ = zeros(Float64, 2, 1)

	α₂ₗ = 1.0 - α₁ₗ
	α₂ᵣ = 1.0 - α₁ᵣ

	j = 1
	(ρᵢₗ[j], Hₜᵢₗ[j], cᵢₗ[j], ∂ρ∂pᵢₗ[j], ∂ρ∂Tᵢₗ[j], ∂Hₜ∂pᵢₗ[j], ∂Hₜ∂Tᵢₗ[j]) = 
	eosNASG(pₗ, uₗ, vₗ, wₗ, Tₗ)
	j = 2
	(ρᵢₗ[j], Hₜᵢₗ[j], cᵢₗ[j], ∂ρ∂pᵢₗ[j], ∂ρ∂Tᵢₗ[j], ∂Hₜ∂pᵢₗ[j], ∂Hₜ∂Tᵢₗ[j]) = 
	eosIdeal(pₗ, uₗ, vₗ, wₗ, Tₗ)
	
	j = 1
	(ρᵢᵣ[j], Hₜᵢᵣ[j], cᵢᵣ[j], ∂ρ∂pᵢᵣ[j], ∂ρ∂Tᵢᵣ[j], ∂Hₜ∂pᵢᵣ[j], ∂Hₜ∂Tᵢᵣ[j]) = 
	eosNASG(pᵣ, uᵣ, vᵣ, wᵣ, Tᵣ)
	j = 2
	(ρᵢᵣ[j], Hₜᵢᵣ[j], cᵢᵣ[j], ∂ρ∂pᵢᵣ[j], ∂ρ∂Tᵢᵣ[j], ∂Hₜ∂pᵢᵣ[j], ∂Hₜ∂Tᵢᵣ[j]) = 
	eosIdeal(pᵣ, uᵣ, vᵣ, wᵣ, Tᵣ)

	ρₗ_ACID = α₁ᵣ * ρᵢₗ[1] + α₂ᵣ * ρᵢₗ[2]
	ρᵣ_ACID = α₁ₗ * ρᵢᵣ[1] + α₂ₗ * ρᵢᵣ[2]

	Y₁ₗ = α₁ᵣ * ρᵢₗ[1] / ρₗ_ACID
	Y₂ₗ = α₂ᵣ * ρᵢₗ[2] / ρₗ_ACID
	Y₁ᵣ = α₁ₗ * ρᵢᵣ[1] / ρₗ_ACID
	Y₂ᵣ = α₂ₗ * ρᵢᵣ[2] / ρₗ_ACID
	
    ∂ρ∂pₗ_ACID = α₁ᵣ * ∂ρ∂pᵢₗ[1] + α₂ᵣ * ∂ρ∂pᵢₗ[2]
    ∂ρ∂pᵣ_ACID = α₁ₗ * ∂ρ∂pᵢᵣ[1] + α₂ₗ * ∂ρ∂pᵢᵣ[2]

	Hₜₗ_ACID = (α₁ᵣ * ρᵢₗ[1] * Hₜᵢₗ[1] + α₂ᵣ * ρᵢₗ[2] * Hₜᵢₗ[2]) / ρₗ_ACID
	Hₜᵣ_ACID = (α₁ₗ * ρᵢᵣ[1] * Hₜᵢᵣ[1] + α₂ₗ * ρᵢᵣ[2] * Hₜᵢᵣ[2]) / ρᵣ_ACID

	∂Hₜ∂pₗ_ACID = ( 
		1.0/ρₗ_ACID*(∂ρ∂pᵢₗ[1]*α₁ᵣ-∂ρ∂pₗ_ACID*Y₁ₗ)*Hₜᵢₗ[1] +
		1.0/ρₗ_ACID*(∂ρ∂pᵢₗ[2]*α₂ᵣ-∂ρ∂pₗ_ACID*Y₂ₗ)*Hₜᵢₗ[2] +
		Y₁ₗ*∂Hₜ∂pᵢₗ[1] + Y₂ₗ*∂Hₜ∂pᵢₗ[2]
	)

	∂Hₜ∂pᵣ_ACID = ( 
		1.0/ρᵣ_ACID*(∂ρ∂pᵢᵣ[1]*α₁ₗ-∂ρ∂pᵣ_ACID*Y₁ᵣ)*Hₜᵢᵣ[1] +
		1.0/ρᵣ_ACID*(∂ρ∂pᵢᵣ[2]*α₂ₗ-∂ρ∂pᵣ_ACID*Y₂ᵣ)*Hₜᵢᵣ[2] +
		Y₁ᵣ*∂Hₜ∂pᵢᵣ[1] + Y₂ᵣ*∂Hₜ∂pᵢᵣ[2]
	)

	∂ρ∂Tₗ_ACID = α₁ᵣ * ∂ρ∂Tᵢₗ[1] + α₂ᵣ * ∂ρ∂Tᵢₗ[2]
	∂ρ∂Tᵣ_ACID = α₁ₗ * ∂ρ∂Tᵢᵣ[1] + α₂ₗ * ∂ρ∂Tᵢᵣ[2]
	
	∂Hₜ∂Tₗ_ACID = ( 
		1.0/ρₗ_ACID*(∂ρ∂Tᵢₗ[1]*α₁ᵣ-∂ρ∂Tₗ_ACID*Y₁ₗ)*Hₜᵢₗ[1] +
		1.0/ρₗ_ACID*(∂ρ∂Tᵢₗ[2]*α₂ᵣ-∂ρ∂Tₗ_ACID*Y₂ₗ)*Hₜᵢₗ[2] +
		Y₁ₗ*∂Hₜ∂Tᵢₗ[1] + Y₂ₗ*∂Hₜ∂Tᵢₗ[2]
	)

	∂Hₜ∂Tᵣ_ACID = ( 
		1.0/ρᵣ_ACID*(∂ρ∂Tᵢᵣ[1]*α₁ₗ-∂ρ∂Tᵣ_ACID*Y₁ᵣ)*Hₜᵢᵣ[1] +
		1.0/ρᵣ_ACID*(∂ρ∂Tᵢᵣ[2]*α₂ₗ-∂ρ∂Tᵣ_ACID*Y₂ᵣ)*Hₜᵢᵣ[2] +
		Y₁ᵣ*∂Hₜ∂Tᵢᵣ[1] + Y₂ᵣ*∂Hₜ∂Tᵢᵣ[2]
	)
	
	∂ρ∂Y₁ₗ_ACID = 0.0#-ρₗ_ACID^2*(1.0/ρᵢₗ[1] + 1.0/ρᵢₗ[2])
	∂ρ∂Y₁ᵣ_ACID = 0.0#-ρᵣ_ACID^2*(1.0/ρᵢᵣ[1] + 1.0/ρᵢᵣ[2])



	return ρₗ_ACID,ρᵣ_ACID,∂ρ∂pₗ_ACID,∂ρ∂pᵣ_ACID, ∂Hₜ∂pₗ_ACID, ∂Hₜ∂pᵣ_ACID,
	Hₜₗ_ACID, Hₜᵣ_ACID, ∂ρ∂Tₗ_ACID, ∂ρ∂Tᵣ_ACID, ∂Hₜ∂Tₗ_ACID, ∂Hₜ∂Tᵣ_ACID,
	∂ρ∂Y₁ₗ_ACID, ∂ρ∂Y₁ᵣ_ACID


end










function eosIdeal(
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	
	cᵥ = 720.0
    γ = 1.4
	
	cₚ = γ*cᵥ
		
	# density of each phase
	ρ = 1.0/( (γ-1.0)*cᵥ*T/p )
	c = √( γ/ρ*p )

	# d(rho)/d(p)
	∂ρ∂p = ρ/p
	
	# d(rho)/d(T)
	∂ρ∂T = -ρ/T

	# d(h)/d(p)
	∂Hₜ∂p = 0.0
	# d(h)/d(T)
	∂Hₜ∂T = γ*cᵥ

	# internal energy of each phase
	# internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	# enthalpy of each phase
	h = γ*cᵥ*T

		   
	# eti = internal_energy + 0.5*usqrt
	Hₜ = h + 0.5*u²

	# cvi = cv
	# cpi = cp	

    return ρ, Hₜ, c, ∂ρ∂p, ∂ρ∂T, ∂Hₜ∂p, ∂Hₜ∂T
end		

function eosNASG(
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	
    p∞ = 621780000.0
	cᵥ = 3610.0
    γ = 1.19
    b = 6.7212e-4
	q = -1177788.0
	
	#= 
    p∞ = 0.0
	cᵥ = 3518.0
    γ = 1.6
    b = 0.0
	q = 0.0
	=#
 
	#=
	# air
    p∞ = 0.0
	cᵥ = 720.0
    γ = 1.4
    b = 0.0
	q = 0.0
	=#
#=
	# helium
    p∞ = 0.0
	cᵥ = 2440.12
    γ = 1.648
    b = 0.0
	q = 0.0
=#
	#=

	# sod water
    p∞ = 6.e8
	cᵥ = 1569.0
    γ = 4.4
    b = 0.0
	q = 0.0

	=#


	cₚ = γ*cᵥ
		
	# density of each phase
	ρ = 1.0/( (γ-1.0)*cᵥ*T/(p+p∞)+b )
	c = √( γ/(ρ*ρ)*(p+p∞)/(1.0/ρ-b) )

	# d(rho)/d(p)
	∂ρ∂p = ρ*ρ*(1.0/ρ-b)/(p+p∞)
	
	# d(rho)/d(T)
	∂ρ∂T = -ρ*ρ*(1.0/ρ-b)/T

	# d(h)/d(p)
	∂Hₜ∂p = b
	# d(h)/d(T)
	∂Hₜ∂T = γ*cᵥ

	# internal energy of each phase
	# internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	# enthalpy of each phase
	h = γ*cᵥ*T + b*p + q

		   
	# eti = internal_energy + 0.5*usqrt
	Hₜ = h + 0.5*u²

	# cvi = cv
	# cpi = cp	

    return ρ, Hₜ, c, ∂ρ∂p, ∂ρ∂T, ∂Hₜ∂p, ∂Hₜ∂T
end