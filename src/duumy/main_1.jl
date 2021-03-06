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
    ð::controls,
    cells::Vector{mesh.Cell}
)

    X = zeros(Float64, length(cells), 1)
    Y = zeros(Float64, length(cells), 8)
    for i in 1:length(cells)
        X[i] = cells[i].x
        Y[i,1] = cells[i].var[ð.p]
        Y[i,2] = cells[i].var[ð.u]
        Y[i,3] = cells[i].var[ð.v]
        Y[i,4] = cells[i].var[ð.T]
        Y[i,5] = cells[i].var[ð.Î±â]
        Y[i,6] = cells[i].var[ð.Ï]
        Y[i,7] = cells[i].var[ð.c]
        Y[i,8] = cells[i].var[ð.Hâ]
        
    end
    plot(X,Y,layout = grid(3, 3), label = ["p" "u" "v" "T" "Î±â" "Ï" "c" "Hâ"] )

    gui()

end


function plotting2D(
    Nx, Ny, 
    ð::controls,
    cells::Vector{mesh.Cell}
)

    #plt = plot(X,VAR,layout = 
    #grid(3, 2),
    #label = ["p" "u" "T" "Yâ" "Ï" "c"] )
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
            VAR1[i,j] = cells[ijk].var[ð.p]
            VAR2[i,j] = cells[ijk].var[ð.Ï]
            #VAR2[i,j] = cells[ijk].var[ð.Î±â]
            VAR3[i,j] = cells[ijk].var[ð.u]
            VAR4[i,j] = cells[ijk].var[ð.v]
            VAR5[i,j] = cells[ijk].var[ð.w]
            VAR6[i,j] = cells[ijk].var[ð.T]

        end
    end

    #plotlyjs()
    #X = 0.5*Îx:Îx:ð.Lx
    #Y = 0.5*Îy:Îy:ð.Ly
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
        x=0.5*Îx:Îx:ð.Lx,#X, # horizontal axis
        y=0.5*Îy:Îy:ð.Ly,#Y, # vertical axis
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
    Ît = 1.e-2
    =#
    
    Nx = 51
    Ny = 1
    Nz = 1
    Lx = 1.0
    Ly = 0.1
    Lz = 0.1
    Ît = 1.e-5

    realMaxIter = 1000000
    pseudoMaxIter = 30
    pseudoMaxResidual = -4.0

    CFL = 0.5
    Lco = 1.0
    Uco = 1.0

    ð = controls(
        Nx,Ny,Nz, Lx,Ly,Lz, 
        realMaxIter,pseudoMaxIter,pseudoMaxResidual, 
        CFL, Ît, Lco, Uco,
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
        ð,
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
    pâ = 0.0
    for cell in cells
        cell.var[ð.p] = 101325.0 - pâ
        cell.var[ð.u] = 0.0
        cell.var[ð.v] = 0.0
        cell.var[ð.w] = 0.0
        cell.var[ð.T] = 300.0
        cell.var[ð.Yâ] = 0.0
        cell.var[ð.Î±â] = 0.0

        if cell.x < 0.4 && cell.y < 0.4
            cell.var[ð.Yâ] = 1.0
            cell.var[ð.Î±â] = 1.0
        end
    end
    =#






    for cell in cells
        cell.var[ð.p] = 101325.0
        cell.var[ð.u] = 1.0
        cell.var[ð.v] = 0.0
        cell.var[ð.w] = 0.0
        cell.var[ð.T] = 300.0
        cell.var[ð.Yâ] = 0.0
        cell.var[ð.Î±â] = 0.0

        if cell.x < 0.1
            cell.var[ð.Yâ] = 1.0
            cell.var[ð.Î±â] = 1.0
        end
    end








    # EOS
    EOS!(ð, cells)
    #EOS_vf!(ð, cells)

    # Transport
    for cell in cells
        cell.var[ð.Î¼] = cell.var[ð.Î±â] * 0.001 + cell.var[ð.Î±â] * 1.e-5
    end

    ð.realIter = 1
    ð.realMaxIter = 100000
    while(
        ð.realIter <= ð.realMaxIter
    )


        # save n-step values
        for cell in cells
            cell.var[ð.pâ¿] = cell.var[ð.p]
            cell.var[ð.uâ¿] = cell.var[ð.u]
            cell.var[ð.vâ¿] = cell.var[ð.v]
            cell.var[ð.wâ¿] = cell.var[ð.w]
            cell.var[ð.Yââ¿] = cell.var[ð.Yâ]
            cell.var[ð.Î±ââ¿] = cell.var[ð.Î±â]
            cell.var[ð.Ïâ¿] = cell.var[ð.Ï]
            cell.var[ð.Hââ¿] = cell.var[ð.Hâ]
        end

        ð.pseudoIter = 1
        ð.pseudoMaxIter = 50
        while(
            ð.pseudoIter â¤ ð.pseudoMaxIter
        )


            totresi,resi1, resi2, resi3 = 
            coupled!(
                ð,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )

            #println(ð.realIter,", ",ð.pseudoIter,", coupled equation success, ",totresi,", ",resi1,", ",resi2,", ",resi3)
            println(ð.realIter,", ",ð.pseudoIter,", coupled equation success, ",totresi)



            resi1 =
            volumefraction_advection!(
                ð,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )
    
            #println(ð.realIter,", ",ð.pseudoIter,", volumefraction_advection success, ",resi1)
    





            # EOS
            EOS!(ð, cells)
            #EOS_vf!(ð, cells)

            # Transport
            for cell in cells
                cell.var[ð.Î¼] = cell.var[ð.Î±â] * 0.001 + cell.var[ð.Î±â] * 1.e-5
            end
            
            # Plotting
            #plotting1D(Nx, Ny, ð, cells)
            #sleep(1.0)


            # Plotting
            plotting1D(Nx, Ny, ð, cells)
            #sleep(1.0)


            ð.pseudoIter += 1

            

        end

        

#=
        resi1 =
        volumefraction!(
            ð,
            cells,
            faces,
            faces_internal,
            faces_boundary,
            faces_boundary_top,
            faces_boundary_bottom,
            faces_boundary_left,
            faces_boundary_right
        )

        println(ð.realIter,", ",ð.pseudoIter,", volumefraction_advection success, ",resi1)


        # Plotting
        plotting1D(Nx, Ny, ð, cells)
        sleep(1.0)
=#

        ð.realIter += 1

    end

    

end






# calculation main
main()




#=
            resi1 =
            momentum!(
                ð,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )

            println(ð.realIter,", ",ð.pseudoIter,", momentum equation success, ",resi1)

            resi1 =
            pressure!(
                ð,
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right
            )

            println(ð.realIter,", ",ð.pseudoIter,", pressure equation success, ",resi1)
        =#

