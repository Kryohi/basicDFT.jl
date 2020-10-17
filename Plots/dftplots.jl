using DataFrames, CSV, Plots

N = 8
rs_Na = 3.93

df8 = DataFrame(CSV.File("./Data/ksfunctions.csv"))

len = count(df8.iteration.==0)
niter = last(df8.iteration)

grid = df8.grid[1:len]
Vext = V_ext.(grid,Rc(8,rs_Na),rho_b(rs_Na))
Vkss = [df8.Vks[i*len+1:(i+1)*len] for i=0:niter]
rhos = [df8.rho[i*len+1:(i+1)*len] for i=0:niter]
eigf_1s = [df8.eigf_1s[i*len+1:(i+1)*len] for i=0:niter]
eigf_1p = [df8.eigf_1p[i*len+1:(i+1)*len] for i=0:niter]

plot(grid,rhos[1:5:end])
plot(grid,Vkss[1:5:end])
