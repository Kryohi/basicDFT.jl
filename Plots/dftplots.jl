using DataFrames, CSV, Plots

rs_Na = 3.93

df8 = DataFrame(CSV.File("./Data/ksfunctions_20.csv"))
en8 = DataFrame(CSV.File("./Data/ksenergy_20.csv"))

len = count(df8.iteration.==0)
ndata = Int(length(df8.iteration)/len)-1
niter = last(df8.iteration)
grid = df8.grid[1:len]
Vext = V_ext.(grid,Rc(8,rs_Na),rho_b(rs_Na))

Vhs = [df8.Vh[i*len+1:(i+1)*len] for i=0:ndata]
Vkss = [df8.Vks[i*len+1:(i+1)*len] for i=0:ndata]
Vxcs = [Vkss[i+1] .- Vhs[i+1] .- Vext for i=0:ndata]
rhos = [df8.rho[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1s = [df8.eigf_1s[i*len+1:(i+1)*len] for i=0:ndata]
eigf_2s = [df8.eigf_2s[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1p = [df8.eigf_1p[i*len+1:(i+1)*len] for i=0:ndata]
eigf_1d = [df8.eigf_1d[i*len+1:(i+1)*len] for i=0:ndata]

plot(grid,rhos[1:1:end])
plot(grid,Vhs[1:1:end])
plot(grid,Vxcs[2:1:end])
plot(grid,Vkss[1:1:end])
plot(grid,eigf_1s[1:1:end])

last_rho = rhos[end]
last_Vks = Vkss[end]
last_Vh = Vhs[end]
last_Vxc = Vxcs[end]
pVh = plot(grid[100:end-100],last_Vh[100:end-100])
savefig(pVh, "./Plots/vh.pdf")

## energy values
plot(0:niter, en8.e1s)
plot(0:niter, en8.E1)
plot(0:niter, en8.E2)
