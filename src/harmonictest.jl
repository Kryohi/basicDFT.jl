include("numerov.jl")

h = 2e-5
grid = 0:h:8-h
V = map(x->0.5*x^2, grid)
nmax = 7

l = 1
bc0 = [0., h^l]

eigv, eigf = Numerov(l, nmax, grid, V, bc_0=bc0, Estep=1e-2)


grid = 0:h:12-h
V = map(x->0.5*(x-6)^2, grid)

eigv, eigf = Numerov(0, nmax, grid, V, Estep=2e-3, verbose=true)
