using DataFrames, CSV, Plots

df = DataFrame(CSV.File("./Data/ksfunctions.csv"))

len = count(df.iteration.==0)
niter = last(df.iteration)

grid = df.grid[1:len]
Vkss = [df.Vks[i*len+1:(i+1)*len] for i=0:niter]
rhos = [df.rho[i*len+1:(i+1)*len] for i=0:niter]
eigf_1s = [df.eigf_1s[i*len+1:(i+1)*len] for i=0:niter]
eigf_1p = [df.eigf_1p[i*len+1:(i+1)*len] for i=0:niter]

plot(grid,rhos[1:2:end])
plot(grid,Vkss[1:2:end])
