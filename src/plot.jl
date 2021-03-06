
residual_global = []
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
        Y[i,5] = cells[i].var[ð.Yâ]
        Y[i,6] = cells[i].var[ð.Ï]
        Y[i,7] = cells[i].var[ð.c]
        Y[i,8] = cells[i].var[ð.Hâ]
        
    end
    l = @layout [grid(2,4) 
                 b{0.2h}  ]
    plt = plot(X,Y,layout = l, 
    legend = false,
    title = ["p" "u" "v" "T" "Yâ" "Ï" "c" "Hâ"],
    label = ["p" "u" "v" "T" "Yâ" "Ï" "c" "Hâ"],
    #formatter = :plain
    digits=2
    )
    push!(residual_global, ð.residual)
    plot!(plt[9], residual_global, 
    legend = false, title = "Residual", label = "Residual")
    plot!(size=(1200,400))

    if length(residual_global) >= 200
        #empty!(residual_global)
        popfirst!(residual_global)
    end


    #p1 = plot(X,Y,layout = grid(4, 2), label = ["p" "u" "v" "T" "Yâ" "Ï" "c" "Hâ"] )
    #plot!(p1[8], residual_global)
    #p2 = plot(ð.time, residual_global)
    #plot(p1, p2)

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
