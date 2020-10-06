#using basicDFT
include("../src/numerov.jl")
using Test, BenchmarkTools

const tol = 1e-4

@testset "basicDFT.jl" begin

    # 1D harmonic oscillator eigenvalue test on Numerov.jl
    nmax = 7
    h = 2e-5    # increase this to speed up the tests, but tol will also have to be increased
    grid = 0:h:10-h
    V = map(x->0.5*(x-5)^2, grid)
    eigv, _ = Numerov(0, nmax, grid, V, Estep=2e-2)
    trueeigv = 0.5:1:6.5
    @test all(x -> abs(x)<tol, eigv .- trueeigv)

    # 3D harmonic oscillator (l=1) eigenvalue tests on Numerov.jl
    V = map(x->0.5*x^2, grid)
    eigv, _ = Numerov(1, nmax, grid, V, bc_0=[0., h^l], Estep=2e-2)
    trueeigv = 2.5:2:14.5
    @test all(x -> abs(x)<tol, eigv .- trueeigv)


    # DFT test on N=8 jellium nanoparticle



end
