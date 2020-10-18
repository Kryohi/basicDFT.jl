using DataFrames, CSV, Plots

rs_Na = 3.93

df8 = DataFrame(CSV.File("./Data/ksfunctions.csv"))
en8 = DataFrame(CSV.File("./Data/ksenergy.csv"))

len = count(df8.iteration.==0)
niter = last(df8.iteration)
grid = df8.grid[1:len]
Vext = V_ext.(grid,Rc(8,rs_Na),rho_b(rs_Na))

Vhs = [df8.Vh[i*len+1:(i+1)*len] for i=0:niter]
Vkss = [df8.Vks[i*len+1:(i+1)*len] for i=0:niter]
Vxcs = [Vkss[i+1] .- Vhs[i+1] .- Vext for i=0:niter]
rhos = [df8.rho[i*len+1:(i+1)*len] for i=0:niter]
eigf_1s = [df8.eigf_1s[i*len+1:(i+1)*len] for i=0:niter]
eigf_2s = [df8.eigf_2s[i*len+1:(i+1)*len] for i=0:niter]
eigf_1p = [df8.eigf_1p[i*len+1:(i+1)*len] for i=0:niter]
eigf_1d = [df8.eigf_1d[i*len+1:(i+1)*len] for i=0:niter]

plot(grid,rhos[1:5:end])
plot(grid,Vhs[1:5:end])
plot(grid,Vxcs[1:5:end])
plot(grid,Vkss[1:5:end])
plot(grid,eigf_2s[1:5:end])


last_Vks = Vkss[end]
plot(grid[800:end],last_Vks[800:end])

## energy values
plot(1:niter, en8.e1s)
plot(1:niter, en8.E1)
plot(1:niter, en8.E2)
