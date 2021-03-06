function oneD_write(
    ð::controls,
    cells::Vector{mesh.Cell},
    save_target
    )


    x = Vector{Float64}()
    y = Array{Float64}(undef, length(cells), 8)
    #X = zeros(Float64, length(cells), 1)
    #Y = zeros(Float64, length(cells), 8)
    for i in 1:length(cells)
        push!(x, cells[i].x)
        y[i,1] = cells[i].var[ð.p]
        y[i,2] = cells[i].var[ð.u]
        y[i,3] = cells[i].var[ð.v]
        y[i,4] = cells[i].var[ð.T]
        y[i,5] = cells[i].var[ð.Yâ]
        y[i,6] = cells[i].var[ð.Ï]
        y[i,7] = cells[i].var[ð.c]
        y[i,8] = cells[i].var[ð.Hâ]
        
    end
    df = DataFrame(x=x, p=y[:,1], u=y[:,2], v=y[:,3], 
    T=y[:,4], Y1=y[:,5], Rho=y[:,6], c=y[:,7], Ht=y[:,8])
    
    CSV.write(save_target, df) 

end


function twoD_write(
    ð::controls,
    cells::Vector{mesh.Cell},
    save_target
    )


    x = Vector{Float64}()
    y = Vector{Float64}()
    val = Array{Float64}(undef, length(cells), 8)
    #X = zeros(Float64, length(cells), 1)
    #Y = zeros(Float64, length(cells), 8)
    for i in 1:length(cells)
        push!(x, cells[i].x)
        push!(y, cells[i].y)
        val[i,1] = cells[i].var[ð.p]
        val[i,2] = cells[i].var[ð.u]
        val[i,3] = cells[i].var[ð.v]
        val[i,4] = cells[i].var[ð.T]
        val[i,5] = cells[i].var[ð.Yâ]
        val[i,6] = cells[i].var[ð.Ï]
        val[i,7] = cells[i].var[ð.c]
        val[i,8] = cells[i].var[ð.Hâ]
        
    end
    df = DataFrame(x=x, y=y, p=val[:,1], u=val[:,2], v=val[:,3], 
    T=val[:,4], Y1=val[:,5], Rho=val[:,6], c=val[:,7], Ht=val[:,8])
    
    CSV.write(save_target, df) 

end