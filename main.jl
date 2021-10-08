using Plots
#using PlotlyJS
using LinearAlgebra
using SparseArrays
using IterativeSolvers
#using AlgebraicMultigrid

using Pardiso

using CSV
using DataFrames

include("./src/structured_grid_uniform.jl")
include("./src/constant.jl")
include("./src/controls.jl")
include("./src/EOS.jl")
include("./src/transport.jl")
include("./src/momentum.jl")
include("./src/pressure.jl")
include("./src/volumefraction.jl")
include("./src/massfraction.jl")
include("./src/energy.jl")
include("./src/coupled_fully.jl")
include("./src/plot.jl")
include("./src/write.jl")
include("./test/oneD_single_phase_subsonic.jl")
include("./test/oneD_single_phase_supsonic.jl")

function main()

    Nx, Ny, Nz, Lx, Ly, Lz, 
    corantNumber, Δt, Δt_steps, Δt_iters, pseudoMaxIter,
    save_time, save_iteration,time_end,
    temporal_discretizationScheme, spatial_discretizationScheme,
    gravity,
    p∞, cᵥ, γ, b, q,
    initial_p, initial_u, initial_v, initial_T, initial_Y,
    left_p_BCtype, left_u_BCtype, left_v_BCtype, left_T_BCtype, left_Y_BCtype,
    right_p_BCtype, right_u_BCtype, right_v_BCtype, right_T_BCtype, right_Y_BCtype,
    bottom_p_BCtype, bottom_u_BCtype, bottom_v_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
    top_p_BCtype, top_u_BCtype, top_v_BCtype, top_T_BCtype, top_Y_BCtype,
    left_p_BCValue, left_u_BCValue, left_v_BCValue, left_T_BCValue, left_Y_BCValue,
    right_p_BCValue, right_u_BCValue, right_v_BCValue, right_T_BCValue, right_Y_BCValue,
    bottom_p_BCValue, bottom_u_BCValue, bottom_v_BCValue, bottom_T_BCValue, bottom_Y_BCValue,
    top_p_BCValue, top_u_BCValue, top_v_BCValue, top_T_BCValue, top_Y_BCValue = 

    Mac_3_test_case()


    realMaxIter = 1000000
    pseudoMaxResidual = -4.0

    CFL = 0.5
    Lco = 1.0
    Uco = 1.0

    👉 = controls(
        Nx,Ny,Nz, Lx,Ly,Lz, 
        temporal_discretizationScheme, spatial_discretizationScheme,
        gravity,
        p∞, cᵥ, γ, b, q,
        left_p_BCtype, left_u_BCtype, left_v_BCtype, left_T_BCtype, left_Y_BCtype,
        right_p_BCtype, right_u_BCtype, right_v_BCtype, right_T_BCtype, right_Y_BCtype,
        bottom_p_BCtype, bottom_u_BCtype, bottom_v_BCtype, bottom_T_BCtype, bottom_Y_BCtype,
        top_p_BCtype, top_u_BCtype, top_v_BCtype, top_T_BCtype, top_Y_BCtype,
        left_p_BCValue, left_u_BCValue, left_v_BCValue, left_T_BCValue, left_Y_BCValue,
        right_p_BCValue, right_u_BCValue, right_v_BCValue, right_T_BCValue, right_Y_BCValue,
        bottom_p_BCValue, bottom_u_BCValue, bottom_v_BCValue, bottom_T_BCValue, bottom_Y_BCValue,
        top_p_BCValue, top_u_BCValue, top_v_BCValue, top_T_BCValue, top_Y_BCValue,
        realMaxIter,pseudoMaxIter,pseudoMaxResidual, 
        corantNumber, CFL, Δt, Lco, Uco,
        0.0, 0, 0, 0.0, 
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
        21,22,23,24,25,26,27,28,29,30,31,
        32,33,34,35,36,37,38,39,40,41
    )
    👉.time = 0.0

    cells = Vector{mesh.Cell}(undef, 0)
    faces = Vector{mesh.Face}(undef, 0)
    faces_internal = Vector{mesh.Face}(undef, 0)
    faces_boundary = Vector{mesh.Face}(undef, 0)
    faces_boundary_top = Vector{mesh.Face}(undef, 0)
    faces_boundary_bottom = Vector{mesh.Face}(undef, 0)
    faces_boundary_left = Vector{mesh.Face}(undef, 0)
    faces_boundary_right = Vector{mesh.Face}(undef, 0)

    structured_grid_uniform!(
        👉,
        cells,
        faces,
        faces_internal,
        faces_boundary,
        faces_boundary_top,
        faces_boundary_bottom,
        faces_boundary_left,
        faces_boundary_right
    )


    # initialization
    for cell in cells
        cell.var[👉.p] = initial_p(cell.x,cell.y)
        cell.var[👉.u] = initial_u(cell.x,cell.y)
        cell.var[👉.v] = initial_v(cell.x,cell.y)
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = initial_T(cell.x,cell.y)
        cell.var[👉.Y₁] = initial_Y(cell.x,cell.y)
        cell.var[👉.α₁] = initial_Y(cell.x,cell.y)

    end

    # EOS
    EOS!(👉, cells)
    #EOS_vf!(👉, cells)

    # Transport
    for cell in cells
        cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
    end
    

    # save n-step values
    for cell in cells
        cell.var[👉.pⁿ] = cell.var[👉.p]
        cell.var[👉.uⁿ] = cell.var[👉.u]
        cell.var[👉.vⁿ] = cell.var[👉.v]
        cell.var[👉.wⁿ] = cell.var[👉.w]
        cell.var[👉.Y₁ⁿ] = cell.var[👉.Y₁]
        cell.var[👉.α₁ⁿ] = cell.var[👉.α₁]
        cell.var[👉.ρⁿ] = cell.var[👉.ρ]
        cell.var[👉.Hₜⁿ] = cell.var[👉.Hₜ]
    end
    

    👉.realIter = 1
    👉.realMaxIter = 1000000
    saveΔt = 👉.Δt
    Δt_iters_save = 1
    while(
        👉.realIter <= 👉.realMaxIter
    )

        println("real-time Step: $(👉.realIter) \t Time: $(👉.time)")


        # save n-1 step values
        for cell in cells
            cell.var[👉.pⁿ⁻¹] = cell.var[👉.pⁿ]
            cell.var[👉.uⁿ⁻¹] = cell.var[👉.uⁿ]
            cell.var[👉.vⁿ⁻¹] = cell.var[👉.vⁿ]
            cell.var[👉.wⁿ⁻¹] = cell.var[👉.wⁿ]
            cell.var[👉.Y₁ⁿ⁻¹] = cell.var[👉.Y₁ⁿ]
            cell.var[👉.α₁ⁿ⁻¹] = cell.var[👉.α₁ⁿ]
            cell.var[👉.ρⁿ⁻¹] = cell.var[👉.ρⁿ]
            cell.var[👉.Hₜⁿ⁻¹] = cell.var[👉.Hₜⁿ]
        end

        # save n-step values
        for cell in cells
            cell.var[👉.pⁿ] = cell.var[👉.p]
            cell.var[👉.uⁿ] = cell.var[👉.u]
            cell.var[👉.vⁿ] = cell.var[👉.v]
            cell.var[👉.wⁿ] = cell.var[👉.w]
            cell.var[👉.Y₁ⁿ] = cell.var[👉.Y₁]
            cell.var[👉.α₁ⁿ] = cell.var[👉.α₁]
            cell.var[👉.ρⁿ] = cell.var[👉.ρ]
            cell.var[👉.Hₜⁿ] = cell.var[👉.Hₜ]
        end

        if Δt_iters[Δt_iters_save] == 👉.realIter
            Δt_iters_save += 1
            👉.Δt = Δt_steps[Δt_iters_save]
        else
            👉.Δt = Δt_steps[Δt_iters_save]
        end
        

        # timestep from corantNumber
        #=
        Δtₘᵢₙ = 1.e9
        for cell in cells
            U = [cell.var[👉.u] cell.var[👉.v] cell.var[👉.w]]
            Δtₘᵢₙ = min(Δtₘᵢₙ, cell.Ω^0.333 / norm(U))
        end
        if 👉.corantNumber*Δtₘᵢₙ > saveΔt
            👉.Δt = 👉.corantNumber*Δtₘᵢₙ
        else
            👉.Δt = saveΔt
        end
        if 👉.realIter == 1
            👉.Δt = saveΔt
        end
        =#


        for ii in 1:5

            resi1 =
            massfraction!(
                👉,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )

            println(👉.realIter,", ",👉.pseudoIter,", volumefraction_advection success, ",resi1)

            # EOS
            EOS!(👉, cells)
            #EOS_vf!(👉, cells)

            # Transport
            for cell in cells
                cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
            end
            

        end




        👉.pseudoIter = 1
        #👉.pseudoMaxIter = 25
        while(
            👉.pseudoIter ≤ 👉.pseudoMaxIter[Δt_iters_save]
        )


            totresi,resi1, resi2, resi3 = 
            coupled!(
                👉,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )

            #println(👉.realIter,", ",👉.pseudoIter,", coupled equation success, ",totresi,", ",resi1,", ",resi2,", ",resi3)
            println(👉.realIter,", ",👉.pseudoIter,", coupled equation success, ",totresi)

            👉.residual = totresi

            # EOS
            EOS!(👉, cells)
            #EOS_vf!(👉, cells)

            # Transport
            for cell in cells
                cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
            end
            
            # Plotting
            if Ny == 1
                plotting1D(Nx, Ny, 👉, cells)
            else
                plotting2D(Nx, Ny, 👉, cells)
            end

            #sleep(1.0)


            👉.pseudoIter += 1

            

        end


        if (👉.time-👉.Δt <= save_time <= 👉.time) ||
            ( 👉.realIter % save_iteration == 0 )
            save_file::String = "save\\" * string(👉.realIter) * 
            "_" * string(round(👉.time; digits=9)) * ".csv"
            if Ny == 1
                oneD_write(👉, cells, save_file)
            else
                twoD_write(👉, cells, save_file)
            end
        end

        if 👉.time > time_end
            break
        end
        
        👉.time += 👉.Δt

        👉.realIter += 1

    end

    

end






# calculation main
main()



