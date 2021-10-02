include("./structured_grid_uniform.jl")
include("./constant.jl")
include("./controls.jl")
include("./EOS.jl")
include("./transport.jl")
include("./momentum.jl")
include("./pressure.jl")
include("./volumefraction.jl")
include("./massfraction.jl")
include("./energy.jl")
include("./coupled_fully.jl")


using Plots
#using PlotlyJS
using LinearAlgebra
using SparseArrays
using IterativeSolvers
#using AlgebraicMultigrid

using Pardiso

function plotting1D(
    Nx, Ny, 
    👉::controls,
    cells::Vector{mesh.Cell}
)

    X = zeros(Float64, length(cells), 1)
    Y = zeros(Float64, length(cells), 8)
    for i in 1:length(cells)
        X[i] = cells[i].x
        Y[i,1] = cells[i].var[👉.p]
        Y[i,2] = cells[i].var[👉.u]
        Y[i,3] = cells[i].var[👉.v]
        Y[i,4] = cells[i].var[👉.T]
        Y[i,5] = cells[i].var[👉.α₁]
        Y[i,6] = cells[i].var[👉.ρ]
        Y[i,7] = cells[i].var[👉.c]
        Y[i,8] = cells[i].var[👉.Hₜ]
        
    end
    plot(X,Y,layout = grid(3, 3), label = ["p" "u" "v" "T" "α₁" "ρ" "c" "Hₜ"] )

    gui()

end


function plotting2D(
    Nx, Ny, 
    👉::controls,
    cells::Vector{mesh.Cell}
)

    #plt = plot(X,VAR,layout = 
    #grid(3, 2),
    #label = ["p" "u" "T" "Y₁" "ρ" "c"] )
    #plot(plt)
    #contourf!(X,Y,VAR)

    X = zeros(Float64, Nx)
    Y = zeros(Float64, Ny)
    VAR1 = zeros(Float64, Nx, Ny)
    VAR2 = zeros(Float64, Nx, Ny)
    VAR3 = zeros(Float64, Nx, Ny)
    VAR4 = zeros(Float64, Nx, Ny)
    VAR5 = zeros(Float64, Nx, Ny)
    VAR6 = zeros(Float64, Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            k=1
            ijk = i + Nx*(j-1) + Nx*Ny*(k-1)
            X[i] = cells[ijk].x
            Y[j] = cells[ijk].y
            VAR1[i,j] = cells[ijk].var[👉.p]
            VAR2[i,j] = cells[ijk].var[👉.ρ]
            #VAR2[i,j] = cells[ijk].var[👉.α₁]
            VAR3[i,j] = cells[ijk].var[👉.u]
            VAR4[i,j] = cells[ijk].var[👉.v]
            VAR5[i,j] = cells[ijk].var[👉.w]
            VAR6[i,j] = cells[ijk].var[👉.T]

        end
    end

    #plotlyjs()
    #X = 0.5*Δx:Δx:👉.Lx
    #Y = 0.5*Δy:Δy:👉.Ly
    #X = repeat(reshape(x, 1, :), length(y), 1)
    #Y = repeat(y, 1, length(x))
    #plot(contour(X, Y, VAR2, fill = true))
    plot(
        heatmap(X, Y, VAR1', c = :bluesreds),
        heatmap(X, Y, VAR2', c = :bluesreds),
        heatmap(X, Y, VAR3', c = :bluesreds),
        heatmap(X, Y, VAR4', c = :bluesreds),
        heatmap(X, Y, VAR5', c = :bluesreds),
        heatmap(X, Y, VAR6', c = :bluesreds);
        layout = grid(3, 2)
    )

    gui()
#=
    plot(contour(
        x=0.5*Δx:Δx:👉.Lx,#X, # horizontal axis
        y=0.5*Δy:Δy:👉.Ly,#Y, # vertical axis
        z=VAR2'#VAR[:,5]'
    ))
=#



end


function main()

    #=
    Nx = 50
    Ny = 50
    Nz = 1
    Lx = 1.0
    Ly = 1.0
    Lz = 0.1
    Δt = 1.e-2
    =#
    
    Nx = 50
    Ny = 50
    Nz = 1
    Lx = 1.0
    Ly = 1.0
    Lz = 0.1
    Δt = 1.e-4

    realMaxIter = 1000000
    pseudoMaxIter = 3
    pseudoMaxResidual = -4.0

    CFL = 0.5
    Lco = 1.0
    Uco = 1.0

    👉 = controls(
        Nx,Ny,Nz, Lx,Ly,Lz, 
        realMaxIter,pseudoMaxIter,pseudoMaxResidual, 
        CFL, Δt, Lco, Uco,
        0.0, 0, 0, 0.0, 
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
        21,22,23,24,25,26,27,28,29,30,31
    )

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

    # dam break
    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 0.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
        cell.var[👉.α₁] = 0.0

        if cell.x < 0.4 && cell.y < 0.4
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end





    #=

    # 1D interface advection with constant velocity
    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 1.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
        cell.var[👉.α₁] = 0.0

        if cell.x < 0.5
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end
    =#

    #=
    # 1D helium-bubble in air
    for cell in cells

        if cell.x < 0.3
            cell.var[👉.p] = 1.245e5
            cell.var[👉.u] = 55.33
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 319.48
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 1.e5
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end

        
        if 0.5 < cell.x < 0.7
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end
=#

#=
    # 1D sod shock
    for cell in cells

        if cell.x < 0.5
            cell.var[👉.p] = 1.0
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.003484
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 0.1
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.002787
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end

    end
=#

#=

    # 1D high-pressure water and low-pressure air shock tube
    for cell in cells

        if cell.x < 0.7
            cell.var[👉.p] = 1.e9
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        else
            cell.var[👉.p] = 1.e5
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 6.968
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end

    end
=#



    # EOS
    EOS!(👉, cells)
    #EOS_vf!(👉, cells)

    # Transport
    for cell in cells
        cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
    end

    👉.realIter = 1
    👉.realMaxIter = 1000000
    while(
        👉.realIter <= 👉.realMaxIter
    )

        if 👉.realIter < 10 
            👉.Δt = 5.e-4
        elseif 👉.realIter < 30
            👉.Δt = 3.e-3
        else
            👉.Δt = 3.e-3
        end
            

    
        println("real-time Step: $(👉.realIter) \t Time: $(👉.time)")

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
            👉.pseudoIter ≤ 👉.pseudoMaxIter
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


            # EOS
            EOS!(👉, cells)
            #EOS_vf!(👉, cells)

            # Transport
            for cell in cells
                cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
            end
            
            # Plotting
            plotting2D(Nx, Ny, 👉, cells)
            #sleep(1.0)


            👉.pseudoIter += 1

            

        end

        
        👉.time += 👉.Δt

        👉.realIter += 1

    end

    

end






# calculation main
main()



