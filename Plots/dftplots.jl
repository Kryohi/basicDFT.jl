using DataFrames, CSV, Plots

nucl = "K"
N = 20

if nucl == "Na"
    rs = 3.93
else
    rs = 4.86
end

df = DataFrame(CSV.File("./Data/ksfunctions_$(nucl)_$N.csv"))
en = DataFrame(CSV.File("./Data/ksenergy_$(nucl)_$N.csv"))

<<<<<<< HEAD

cos_single(x,l,c) = cos((x-c)*pi/l)^2 * (abs((x-c)*pi/l)<pi/2) + 1e-9
rho = cos_single.(grid, 4, 4.9)
rho = rho .* 8 ./ norm(rho,1)
plot(grid,rho)

Vhs = [df8.Vh[i*len+1:(i+1)*len] for i=0:ndata]
Vkss = [df8.Vks[i*len+1:(i+1)*len] for i=0:ndata]
=======
len = count(df.iteration.==0)
ndata = Int(length(df.iteration)/len)-1
niter = last(df.iteration)
grid = df.grid[1:len]
Vext = V_ext.(grid, Rc(N,rs), rho_b(rs))

Vhs = [df.Vh[i*len+1:(i+1)*len] for i=0:ndata]
Vkss = [df.Vks[i*len+1:(i+1)*len] for i=0:ndata]
>>>>>>> 9e4aa5b714f26c468117d8710a33c65c965a5595
Vxcs = [Vkss[i+1] .- Vhs[i+1] .- Vext for i=0:ndata]
rhos = [df.rho[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1s = [df.eigf_1s[i*len+1:(i+1)*len] for i=0:ndata]
eigf_2s = [df.eigf_2s[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1p = [df.eigf_1p[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1d = [df.eigf_1d[i*len+1:(i+1)*len] for i=0:ndata]

plot(grid,rhos[1:1:end])
plot(grid,Vhs[1:1:end], ylims=(0, 0.8))
plot(grid,Vxcs[1:1:end])
plot(grid,Vkss[1:1:end], ylims=(-Inf, 0.8))
plot(grid,eigf_2s[1:1:end])

# plot of delta


## Plots of the final iteration

last_rho = rhos[end]
last_Vh = Vhs[end]
last_Vxc = Vxcs[end]
last_Vks = Vkss[end]
last_1s = eigf_1s[end]
last_1p = eigf_1p[end]
last_1d = eigf_1d[end]
last_2s = eigf_2s[end]

# plots of the potential components
plot(grid, Vext, legend=false, xlabel="r", ylabel="Vext")
savefig("./Plots/Vext_$(nucl)_$N.pdf")
plot(grid[100:end-100], last_Vh[100:end-100], legend=false)
savefig("./Plots/Vh_$(nucl)_$N.pdf")
plot(grid, last_Vxc, legend=false)
savefig("./Plots/Vxc_$(nucl)_$N.pdf")
plot(grid[100:end-100], last_Vks[100:end-100], legend=false)
savefig("./Plots/Vks_$(nucl)_$N.pdf")

# numerical instability
plot(grid[29900:30000],last_Vh[29900:30000])


# plot of the radial wavefunctions
plot(grid[1:end÷2], last_1s[1:end÷2], label="1s")
plot!(grid[1:end÷2], last_1p[1:end÷2], label="1p")
if N==20
    plot!(grid[1:end÷2], last_1d[1:end÷2], label="1d")
    plot!(grid[1:end÷2], last_2s[1:end÷2], label="2s")
end
savefig("./Plots/eigf_$(nucl)_$N.pdf")

# plot of the radial probabilities
plot(grid[1:end÷2], last_1s[1:end÷2].^2, label="1s")
plot!(grid[1:end÷2], last_1p[1:end÷2].^2, label="1p")
if N==20
    plot!(grid[1:end÷2], last_1d[1:end÷2].^2, label="1d")
    plot!(grid[1:end÷2], last_2s[1:end÷2].^2, label="2s")
end
savefig("./Plots/eigf_squared_$(nucl)_$N.pdf")

# plot of the total electron radial density
plot(grid[1:end÷2], last_rho[1:end÷2])
savefig("./Plots/rho_$(nucl)_$N.pdf")



## energy values
plot(0:niter, en.e1s)
plot(0:niter, en.E1)
plot(0:niter, en.E2)
