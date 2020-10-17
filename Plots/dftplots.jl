using DataFrames, CSV, Plots

df = DataFrame(CSV.File("./Data/ksfunctions.csv"))

len = count(df.iteration.==0)
niter = last(df.iteration)

grid = df.grid[1:len]
Vkss = [df.Vks[i*len+1:(i+1)*len] for i=0:niter]
rhos = [df.rho[i*len+1:(i+1)*len] for i=0:niter]
eigf_1s = [df.eigf_1s[i*len+1:(i+1)*len] for i=0:niter]
eigf_1p = [df.eigf_1p[i*len+1:(i+1)*len] for i=0:niter]

plot(grid,rhos[1:3:end])
plot(grid,Vkss[1:3:end])





cos_single(x,l,c) = cos((x-c)*pi/l)^2 * (abs((x-c)*pi/l)<pi/2) + 1e-9
rho = cos_single.(grid, 4, 4.9)
rho = rho .* 8 ./ norm(rho,1)
plot(grid,rho)
