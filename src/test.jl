function construct_A_matrix_explicit!(
    ๐::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell}
)

    for i in 1:length(cells)

        ฯ = cells[i].var[๐.ฯ]
        u = cells[i].var[๐.u]
        v = cells[i].var[๐.v]
        w = cells[i].var[๐.w]
        T = cells[i].var[๐.T]
        Yโ = cells[i].var[๐.Yโ]
        Hโ = cells[i].var[๐.Hโ]
        โฯโp = cells[i].var[๐.โฯโp]
        โฯโT = cells[i].var[๐.โฯโT]
        โHโโp = cells[i].var[๐.โHโโp]
        โHโโT = cells[i].var[๐.โHโโT]
        โฯโYโ = cells[i].var[๐.โฯโYโ]
        โHโโYโ = cells[i].var[๐.โHโโYโ]

        T = zeros(Float64, 6, 6)
        T[1, 1] = โฯโp
        T[2, 1] = โฯโp*u
        T[3, 1] = โฯโp*v
        T[4, 1] = โฯโp*w
        T[5, 1] = โฯโp*Hโ + ฯ*โHโโp - 1.0
        
        T[2, 2] = ฯ
        T[5, 2] = ฯ*u
        
        T[3, 3] = ฯ
        T[5, 3] = ฯ*v
        
        T[4, 4] = ฯ
        T[5, 4] = ฯ*w
        
        T[1, 5] = โฯโT
        T[2, 5] = โฯโT*u
        T[3, 5] = โฯโT*v
        T[4, 5] = โฯโT*w
        T[5, 5] = โฯโT*Hโ + ฯ*โHโโT

        T[6, 1] = โฯโp*Yโ
        T[6, 5] = โฯโT*Yโ

        T[6, 6] = โฯโYโ*Yโ
        T[6, 6] += ฯ

        T[1, 6] = โฯโYโ
        T[2, 6] = โฯโYโ*u
        T[3, 6] = โฯโYโ*v
        T[4, 6] = โฯโYโ*w
        T[5, 6] = โฯโYโ*Hโ + ฯ*โHโโYโ

        P = deepcopy(T)
        ฮฒ = 1.0/cells[i].var[๐.Vแตฃ]^2.0
            + โฯโp - 1.0/cells[i].var[๐.c]^2.0
        P[1, 1] = ฮฒ
        P[2, 1] = ฮฒ*u
        P[3, 1] = ฮฒ*v
        P[4, 1] = ฮฒ*w
        P[5, 1] = ฮฒ*Hโ + ฯ*โHโโp - 1.0
        P[6, 1] = ฮฒ*Yโ

        A[i, :, :] = T[:, :]./๐.ฮt


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
    ๐::controls, 
    RHS::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
)


    for face in faces_internal
        cpi2 = 0.0
		wโ = 0.0
		pโ = face.varโ[๐.p]
		ฯโ = face.varโ[๐.ฯ]
		cโ = face.varโ[๐.c]
		pแตฃ = face.varแตฃ[๐.p]
		ฯแตฃ = face.varแตฃ[๐.ฯ]
		cแตฃ = face.varแตฃ[๐.c]
		preLs = abs(pโ) + 0.1 * ฯโ*cโ*cโ
		preRs = abs(pแตฃ) + 0.1 * ฯแตฃ*cแตฃ*cแตฃ
		wโ = min(preLs/preRs,preRs/preLs)

        Uโ = 0.0

        flux = zeros(Float64, 6, 1)
        flux = KT_KNP(
            face.varโ[๐.p],face.varโ[๐.u],face.varโ[๐.v],face.varโ[๐.w],
            face.varโ[๐.T],face.varโ[๐.Yโ],face.varโ[๐.ฯ],face.varโ[๐.Hโ],face.varโ[๐.c],
            face.varแตฃ[๐.p],face.varแตฃ[๐.u],face.varแตฃ[๐.v],face.varแตฃ[๐.w],
            face.varแตฃ[๐.T],face.varแตฃ[๐.Yโ],face.varแตฃ[๐.ฯ],face.varแตฃ[๐.Hโ],face.varแตฃ[๐.c],
            ๐.Lco,๐.Uco,๐.ฮt,
            wโ, wโ, cpi2, Uโ,
            face.nฬ[1],face.nฬ[2],face.nฬ[3]
            )


        RHS[face.owner, :] -= flux[:]*face.ฮS/cells[face.owner].ฮฉ #* ๐.ฮt
        RHS[face.neighbour, :] += flux[:]*face.ฮS/cells[face.neighbour].ฮฉ #* ๐.ฮt

    end


    for face in faces_boundary

		Uโ = 0.0

        flux = zeros(Float64, 6, 1)
        flux = KT_KNP(
            face.varโ[๐.p],face.varโ[๐.u],face.varโ[๐.v],face.varโ[๐.w],
            face.varโ[๐.T],face.varโ[๐.Yโ],face.varโ[๐.ฯ],face.varโ[๐.Hโ],face.varโ[๐.c],
            face.varแตฃ[๐.p],face.varแตฃ[๐.u],face.varแตฃ[๐.v],face.varแตฃ[๐.w],
            face.varแตฃ[๐.T],face.varแตฃ[๐.Yโ],face.varแตฃ[๐.ฯ],face.varแตฃ[๐.Hโ],face.varแตฃ[๐.c],
            ๐.Lco,๐.Uco,๐.ฮt,
            1.0, 1.0, 1.0, Uโ,
            face.nฬ[1],face.nฬ[2],face.nฬ[3]
            )

        RHS[face.owner, :] -= flux[:]*face.ฮS/cells[face.owner].ฮฉ #* ๐.ฮt

    end

	

end






function KT_KNP(
    pโ,uโ,vโ,wโ,Tโ,Yโโ,ฯโ,Hโโ,cโ,
    pแตฃ,uแตฃ,vแตฃ,wแตฃ,Tแตฃ,Yโแตฃ,ฯแตฃ,Hโแตฃ,cแตฃ,
    Lco,Uco,ฮt,
    wโ,wโ,cpi2, Uโ,
    nx,ny,nz
)

	Uโโ = uโ*nx + vโ*ny + wโ*nz
	Uโแตฃ = uแตฃ*nx + vแตฃ*ny + wแตฃ*nz

	Fmax =  max(max(Uโโ+cโ,Uโแตฃ+cแตฃ),0.0)
	Fmin = -min(min(Uโโ-cโ,Uโแตฃ-cแตฃ),0.0)

	ฮฑโบ = Fmax / (Fmax+Fmin)
	ฮฑโป = Fmin / (Fmax+Fmin)
	ฮฑโบโป = Fmax*Fmin / (Fmax+Fmin)

	ฯโบ = ฯโ * (ฮฑโบ * Uโโ + ฮฑโบโป)
	ฯโป = ฯแตฃ * (ฮฑโป * Uโแตฃ - ฮฑโบโป)

	mฬโ = ฯโบ
	mฬแตฃ = ฯโป
	
	#mฬ = ฯโบ + ฯโป
	#mฬโ = 0.5*(mฬ+abs(mฬ))
	#mฬแตฃ = 0.5*(mฬ-abs(mฬ))

	
	pโแตฃ = 0.5*(pโ + pแตฃ)

	# comp. convective flux
	flux = zeros(Float64,6,1)
	flux[1] = mฬโ + mฬแตฃ
	flux[2] = mฬโ*uโ + mฬแตฃ*uแตฃ + pโแตฃ*nx
	flux[3] = mฬโ*vโ + mฬแตฃ*vแตฃ + pโแตฃ*ny
	flux[4] = mฬโ*wโ + mฬแตฃ*wแตฃ + pโแตฃ*nz
	flux[5] = mฬโ*Hโโ + mฬแตฃ*Hโแตฃ
	flux[6] = mฬโ*Yโโ + mฬแตฃ*Yโแตฃ

	return flux


end








