using DataFrames, CSV, Plots, ColorSchemes

Rc(N,rs) = cbrt(N)*rs # radius of the positive jellium
rho_b(N,rs) = 3*N/(4π*rs^3) # density of charge inside the nucleus

nucl = "Na"
N = 20

if nucl == "Na"
    rs = 3.93
else
    rs = 4.86
end

df = DataFrame(CSV.File("./Data/ksfunctions_$(nucl)_$N.csv"))
en = DataFrame(CSV.File("./Data/ksenergy_$(nucl)_$N.csv"))

@show len = count(df.iteration.==0)
@show ndata = Int(length(df.iteration)/len)-1
@show niter = last(df.iteration)
s = 1 # we take a point every s, to have more lightweight pdf files
grid = df.grid[1:s:len]
Vext = V_ext.(grid, Rc(N,rs), rho_b(N,rs))

Vhs = [df.Vh[i*len+1:s:(i+1)*len] for i=0:ndata]
Vkss = [df.Vks[i*len+1:s:(i+1)*len] for i=0:ndata]
Vxcs = [Vkss[i+1] .- Vhs[i+1] .- Vext for i=0:ndata]
rhos = [df.rho[i*len+1:s:(i+1)*len] for i=0:ndata]
eigf_1s = [df.eigf_1s[i*len+1:s:(i+1)*len] for i=0:ndata]
eigf_2s = [df.eigf_2s[i*len+1:s:(i+1)*len] for i=0:ndata]
eigf_1p = [df.eigf_1p[i*len+1:s:(i+1)*len] for i=0:ndata]
eigf_1d = [df.eigf_1d[i*len+1:s:(i+1)*len] for i=0:ndata]


## plots of the evolution of the functions
N_plots = 10
strd = 10#Int(floor(ndata/N_plots))
col = palette([:orange, :purple], length(rhos[1:strd:end]))
xlim = floor(Int,grid[end]*1)

plot(grid,rhos[1:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_rho_$(nucl)_$N.pdf")
plot(grid,Vhs[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_Vh_$(nucl)_$N.pdf")
plot(grid,Vxcs[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_Vxc_$(nucl)_$N.pdf")
plot(grid,Vkss[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_Vks_$(nucl)_$N.pdf")
plot(grid,eigf_1s[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_1s_$(nucl)_$N.pdf")
plot(grid,eigf_1p[1:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_1p_$(nucl)_$N.pdf")
plot(grid,eigf_1d[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_1d_$(nucl)_$N.pdf")
plot(grid,eigf_2s[2:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_2s_$(nucl)_$N.pdf")



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
plot(grid[1:end], Vext[1:end], legend=false, xlabel="r", ylabel="Vext")
savefig("./Plots/Vext_$(nucl)_$N.pdf")
plot(grid, last_Vh, legend=false)
savefig("./Plots/Vh_$(nucl)_$N.pdf")
plot(grid[1:end], last_Vxc[1:end], legend=false)
savefig("./Plots/Vxc_$(nucl)_$N.pdf")
plot(grid[1:end-500], last_Vks[1:end-500], legend=false)
savefig("./Plots/Vks_$(nucl)_$N.pdf")

# numerical instability
#plot(grid[29900:30000],last_Vh[29900:30000])
plot(grid[17500:17900], last_Vks[17500:17900], legend=false)


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
plot(en.e1s, label="1s")
plot!(en.e1p, label="1p")
plot!(en.e1d, label="1d")
plot!(en.e2s, label="2s")
plot(en.E1)
plot(en.E2)


df = DataFrame(CSV.File("./Data/Vh.csv"))
Vhq = df.Vh
Vksq = df.Vks

plot(grid, Vhq)
plot(grid[1000:end-910], Vksq[1000:end-910])

df = DataFrame(CSV.File("./Data/rho.csv"))
plot(grid, df.rho)
plot(grid, df.eig1s)
plot(grid, df.eig2s)

# initial conditions
rmax = 32
h = 1e-4
N = 20
grid = Vector(h:h:rmax)
cos_single(x,w,c) = cos((x-c)*π/w)^4 * (abs((x-c)*π/w) < π/2) + 1e-4 + (1e-4*(x≥c+w/2)*exp(-x+c-w/2) - 1e-4) + 1e-18
smoothed_theta(x,x1,x2,w) =  (x>x1 && x<x2) + (x≥x2)*exp(-(x-x2)/w) + (x≤x1)*exp((x-x1)/w)
rho = cos_single.(grid, 100, 20)
#rho = 8 ./(20 .+ (grid .- 20).^4)
#rho = smoothed_theta.(grid,0,20,1)
rho = last_rho
rho = rho .* N ./ simpson_integral(rho, h) #norm(rho,1)#
Vext = V_ext.(grid, Rc(N,rs_Na), rho_b(N,rs_Na))
Vh = 1 .* V_h(grid, rho) ./100
Vxc = V_xc(rho.*grid.^2)
Vks = Vext .+ Vh .+ Vxc
bc_0 = [rho[1:2]./N; rho[1:2]./N; rho[1:2]./N].^(0.5)
bc_end = [exp(-grid[end]); exp(-grid[end-1]); exp(-grid[end]); exp(-grid[end-1]); exp(-grid[end]); exp(-grid[end-1])]
# calculation of the eigenfunctions with the effective mean-field potential
eigv_l0, eigf_l0 = Numerov(0, 2, grid, Vks, bc_0=bc_0[1:2], bc_end=bc_end[1:2],
   Estep=1e-2, maxiter=10^4, verbose=false)
eigv_l1, eigf_l1 = Numerov(1, 1, grid, Vks, bc_0=bc_0[3:4], bc_end=bc_end[3:4],
   Estep=1e-2, maxiter=10^4, verbose=false)
eigv_l2, eigf_l2 = Numerov(2, 1, grid, Vks, bc_0=bc_0[5:6], bc_end=bc_end[5:6],
   Estep=1e-2, maxiter=10^4, verbose=true)
rho = 2 .* eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1].^2
if N>8
    rho .+= 2 .* eigf_l0[:,2].^2 + 10 .* eigf_l2[:,1].^2
end

plot(grid, Vxc)
plot(grid, Vext)
plot(grid, Vh)
plot(grid, Vks)
plot(grid, rho)

plot(grid[1:end÷2], eigf_l0[1:end÷2,1], label="1s")
plot!(grid[1:end÷2], eigf_l1[1:end÷2,1], label="1p")
if N==20
    plot!(grid[1:end÷2], eigf_l2[1:end÷2,1], label="1d")
    plot!(grid[1:end÷2], eigf_l0[1:end÷2,2], label="2s")
end
