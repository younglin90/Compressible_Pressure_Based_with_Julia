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
include("./coupled_flows.jl")


using Plots
#using PlotlyJS
using LinearAlgebra
using SparseArrays
using IterativeSolvers

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
    
    Nx = 51
    Ny = 1
    Nz = 1
    Lx = 1.0
    Ly = 0.1
    Lz = 0.1
    Δt = 1.e-5

    realMaxIter = 1000000
    pseudoMaxIter = 30
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

    
    #=
    p∞ = 0.0
    for cell in cells
        cell.var[👉.p] = 101325.0 - p∞
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
    =#






    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 1.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
        cell.var[👉.α₁] = 0.0

        if cell.x < 0.1
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end








    # EOS
    EOS!(👉, cells)
    #EOS_vf!(👉, cells)

    # Transport
    for cell in cells
        cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
    end

    👉.realIter = 1
    👉.realMaxIter = 100000
    while(
        👉.realIter <= 👉.realMaxIter
    )


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

        👉.pseudoIter = 1
        👉.pseudoMaxIter = 50
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



            resi1 =
            volumefraction_advection!(
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
    
            #println(👉.realIter,", ",👉.pseudoIter,", volumefraction_advection success, ",resi1)
    





            # EOS
            EOS!(👉, cells)
            #EOS_vf!(👉, cells)

            # Transport
            for cell in cells
                cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5
            end
            
            # Plotting
            #plotting1D(Nx, Ny, 👉, cells)
            #sleep(1.0)


            # Plotting
            plotting1D(Nx, Ny, 👉, cells)
            #sleep(1.0)


            👉.pseudoIter += 1

            

        end

        

#=
        resi1 =
        volumefraction!(
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


        # Plotting
        plotting1D(Nx, Ny, 👉, cells)
        sleep(1.0)
=#

        👉.realIter += 1

    end

    

end






# calculation main
main()




#=
            resi1 =
            momentum!(
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

            println(👉.realIter,", ",👉.pseudoIter,", momentum equation success, ",resi1)

            resi1 =
            pressure!(
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

            println(👉.realIter,", ",👉.pseudoIter,", pressure equation success, ",resi1)
        =#

