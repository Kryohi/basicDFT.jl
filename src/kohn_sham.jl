using LinearAlgebra, Printf, DataFrames, CSV
include("functionals.jl")
include("numerov.jl")

# TODO
# fix energy calculation
# compare results with different initial conditions
# check numerical instability of Vh
# save partial data to csv


function solve_KS(N, α, grid, Vext; max_iter=80, stride=1, verbose=false)

      # set the initial trial electron density
      cos_single(x,l,c) = cos((x-c)*π/l)^2 * (abs((x-c)*π/l) < π/2) + 1e-12
      rho = cos_single.(grid, 26, 2)
      rho = rho .* N ./ simpson_integral(rho, h) #norm(rho,1)#
      @show simpson_integral(rho, h)

      # initial boundary conditions for the wavefunctions (s,p,d)
      # the bc for different values of l will become different, so here we concatenate them
      bc_0 = [rho[1:2]./N; rho[1:2]./N; rho[1:2]./N].^(0.5)
      bc_end = [exp(-grid[end]); exp(-grid[end-1]); exp(-grid[end]); exp(-grid[end-1]); exp(-grid[end]); exp(-grid[end-1])]

      # initialize the dataframes to save the data
      data = DataFrame(iteration = zeros(Int16,length(grid)),
                  grid = grid,
                  Vh = zeros(Float64,length(grid)),
                  Vks = zeros(Float64,length(grid)),
                  rho = rho,
                  eigf_1s = zeros(Float64,length(grid)),
                  eigf_2s = zeros(Float64,length(grid)),
                  eigf_1p = zeros(Float64,length(grid)),
                  eigf_1d = zeros(Float64,length(grid)))

      energies = DataFrame(e1s = Float64[],
                        e2s = Float64[],
                        e1p = Float64[],
                        e1d = Float64[],
                        E1 = Float64[],
                        E2 = Float64[])

      @printf("\n\nStarting KS algorithm with N = %d, α = %0.2f\n", N, α)

      # Start of the self-consistent Kohn-Sham method
      for t = 1:max_iter

            @printf("\nITERATION %d\n", t)

            data_step, energy_step = kohn_sham_step(grid, Vext, rho, bc_0, bc_end, N=N, verbose=verbose)

            @show simpson_integral(data_step.rho, h)

            # save partial results to data
            replace!(data_step.iteration, -1 => Int16(t))
            ((t%stride==1) || stride == 1) && append!(data,data_step)
            push!(energies, energy_step)

            # next boundary conditions will be based on the current eigenfunctions
            # NOTE conditions at the end can be either these or a simple decaying exponential, they should both work
            bc_0_new = [data_step.eigf_1s[1:2,1]; data_step.eigf_1p[1:2,1]; data_step.eigf_1d[1:2,1]]
            bc_end_new = [data_step.eigf_1s[end-1:end,1]; data_step.eigf_1p[end-1:end,1]; data_step.eigf_1d[end-1:end,1]]

            # mixing of the eigenfunction's extremes with the old ones
            @show bc_0 .=  α.*bc_0_new .+ (1-α).*bc_0
            @show bc_end .=  α.*bc_end_new .+ (1-α).*bc_end
            @show bc_0 .= bc_0_new
            @show bc_end .= bc_end_new
            #bc_end .=  [-1.;-1.;-1.;-1.;-1.;-1.] # use default exponential in Numerov

            # Check on the convergence by looking at how different is the new density
            @show delta = simpson_integral((data_step.rho .- rho).^2, h)

            # save the found density function to be used and compared in the next iteration
            rho = data_step.rho

            if delta < 1e-4
                  @printf("\nConvergence reached after %d steps with δ = %0.9f\n", t, delta)
                  break
            end
            if t == max_iter
                  @warn("\nExiting from the loop after $max_iter iterations, delta = $delta\n")
            end
      end

      return data, energies
end


# computes the mean-field potential Vks, solves the Shrodinger equation for the
# relevant quantum numbers, computes the new electronic density and returns everything
# as a dataframe
function kohn_sham_step(grid::Vector, Vext::Vector, rho::Vector, bc_0::Vector, bc_end::Vector; N=8, Vks_cutoff=1e5, Estep=1e-3, verbose=false)

      # Hartree potential term, hotspot of the code
      Vh = 1 .* V_h(grid, rho) ./800
      @show minimum(Vh)
      # Kohn-Sham potential, function of rho
      Vks = Vext .+ Vh .+ V_xc(rho)
      # sharp cutoff on the potential, does not have any effect,
      # only used to get results even in the case of a wrong potential calculation
      Vks[Vks.>Vks_cutoff] .= Vks_cutoff
      #Vks .+= Vwall(grid,Vks)

      df = DataFrame(Vh = Vh, Vks = Vks)
      CSV.write("./Data/Vh.csv", df)

      # calculation of the eigenfunctions with the effective mean-field potential
      eigv_l0, eigf_l0 = Numerov(0, 2, grid, Vks, bc_0=bc_0[1:2], bc_end=bc_end[1:2],
       Estep=Estep, maxiter=10^4, verbose=verbose)
      eigv_l1, eigf_l1 = Numerov(1, 1, grid, Vks, bc_0=bc_0[3:4], bc_end=bc_end[3:4],
       Estep=Estep, maxiter=10^4, verbose=verbose)
      eigv_l2, eigf_l2 = Numerov(2, 1, grid, Vks, bc_0=bc_0[5:6], bc_end=bc_end[5:6],
       Estep=Estep, maxiter=10^4, verbose=verbose)

      @show eigv_l0[1]
      @show eigv_l1[1]
      #@show simpson_integral(eigf_l0[:,1].^2, h)
      #@show simpson_integral(eigf_l1[:,1].^2, h)

      # compute the total electron density
      rho_new = 2 .* eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1].^2 .+ 1e-18
      if N>8
            rho_new .+= 2 .* eigf_l0[:,2].^2 + 10 .* eigf_l2[:,1].^2
      end
      @show simpson_integral(rho_new, h)
      rho_new = rho_new .* N ./ simpson_integral(rho_new, h)

      # mixing of the new density with the old one
      rho_new =  α.*rho_new .+ (1-α).*rho

      rho_new = rho_new .* N ./ simpson_integral(rho_new, h)

      df = DataFrame(rho = rho_new, eig1s = eigf_l0[:,1])
      CSV.write("./Data/rho.csv", df)

      # save the computed functions
      data_tmp = DataFrame(iteration = -1 .* ones(Int16,length(grid)),
                        grid = grid,
                        Vh = Vh,
                        Vks = Vks,
                        rho = rho_new,
                        eigf_1s = eigf_l0[:,1],
                        eigf_2s = eigf_l0[:,2],
                        eigf_1p = eigf_l1[:,1],
                        eigf_1d = eigf_l2[:,1])

      # Consistency check through the energy
      @show T_S(grid,rho)
      @show E_ext(grid,rho,Vext)
      @show E_H(grid,rho,Vh)
      @show E_XC(grid,rho)
      @show E1 = T_S(grid,rho) + E_ext(grid,rho,Vext) + E_H(grid,rho,Vh) + E_XC(grid,rho)
      # sum of the eigenvalues - 1/2 hartree energy - exchange
      E2 = eigv_l0[1]*2 + eigv_l1[1]*6 - E_H(grid,rho,Vh) - E_XC(grid,rho)
      if N>8
            E2 += eigv_l0[2]*2 + eigv_l2[1]*10
      end

      verbose && @printf("E_ks = %f\tE_eig = %f\tdiff = %f\n", E1, E2, E2-E1)

      return data_tmp, [eigv_l0; eigv_l1[1]; eigv_l2[1]; E1; E2]
end

function Vwall(grid,Vks)
      Vw = zeros(length(Vks))
      Vw[end-999:end] = (grid[1:1000].*5).^4 .+ Vks[end-999]

      return Vw
end
