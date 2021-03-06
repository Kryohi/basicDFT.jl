using DataFrames, CSV, Plots, ColorSchemes

α = 0.01 # mixing coefficient of the functions
Na = 10.0

df = DataFrame(CSV.File("./Data/gpfunctions_$(α)_$Na.csv"))
en = DataFrame(CSV.File("./Data/gpenergy_$(α)_$Na.csv"))

@show len = count(df.iteration.==0)
@show ndata = Int(length(df.iteration)/len)-1
@show niter = last(df.iteration)
s = 10 # we take a point every s, to have more lightweight pdf files
grid = df.grid[1:s:len]
Vext =  0.5 .* (grid.-0).^2

Vints = [df.Vint[i*len+1:s:(i+1)*len] for i=0:ndata]
phis = [df.phi[i*len+1:s:(i+1)*len] for i=0:ndata]


## plots of the evolution of the functions
N_plots = 20
strd = Int(floor(ndata/N_plots))
col = palette([:orange, :purple], length(phis[1:strd:end]))
xlim = ceil(Int,grid[end]*2/3)

plot(grid, Vints[1:strd:end], ylims=(0,5), palette = col, legend=false)
savefig("./Plots/evol_Vint_$(α)_$Na.pdf")
plot(grid, phis[1:strd:end], xlims=(0,xlim), palette = col, legend=false)
savefig("./Plots/evol_phi_$(α)_$Na.pdf")



# plot of delta


## Plots of the final iteration

last_phi = phis[end]
last_Vint = Vints[end]

# plots of the potential components
plot(grid[1:end], Vext[1:end], legend=false, xlabel="r", ylabel="Vext")
savefig("./Plots/Vext_$(α)_$Na.pdf")
plot(grid[2:6000], last_Vint[2:6000], legend=false)
savefig("./Plots/Vint_$(α)_$Na.pdf")
plot(grid[2:end], last_Vint[2:end].+Vext[2:end], ylims=(0,40), legend=false)
savefig("./Plots/Vtot_$(α)_$Na.pdf")

# numerical instability
#plot(grid[29900:30000],last_Vh[29900:30000])


# plot of the radial wavefunctions
plot(grid[1:end], last_phi[1:end], label="1s")
savefig("./Plots/phi_$(α)_$Na.pdf")

# plot of the radial probabilities
plot(grid[1:end], last_phi[1:end].^2, label="1s")
savefig("./Plots/phi_squared_$(α)_$Na.pdf")


## energy values
plot(1:niter-2, en.e1s)
plot(1:niter-1, en.E1)
plot(1:niter-1, en.E2)

plot(en.e1s[1:end])


## varying α


cos_single(x,l,c) = cos((x-c)*π/l)^2 * (abs((x-c)*π/l) < π/2) + 1e-12
phi = cos_single.(grid, 8, -1)
phi = phi ./ sqrt(sum(phi.^2 .* h))

Na = 1000
Vint(x,phi) = 4π*Na*(phi^2)/x^2# .+ 1e-11))
Vext = 0.5 .* (grid .- 0).^2

Vint.(grid,phi)
4π*Na.*phi.^2 ./ grid.^2
plot(grid[1:10000], phi[1:10000].^2)
plot(grid[1:1000], Vint.(grid,phi)[1:1000])

norm2 = sum(phi.^2 .* h)
