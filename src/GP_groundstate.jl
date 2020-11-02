using LinearAlgebra, Printf, DataFrames, CSV
include("functionals.jl")
include("numerov.jl")


function solve_GP(Na, α, grid::Vector, Vext::Vector; Estep=2e-3, max_iter=500, stride=1, verbose=false)

      # set the initial trial wavefunction
      cos_single(x,l,c) = cos((x-c)*π/l)^2 * (abs((x-c)*π/l) < π/2) + 1e-15
      phi = cos_single.(grid, 8, 3.)
      phi = phi ./ norm(phi)

      # initial boundary condition for the wavefunction
      bc_0 = phi[1:2]
      bc_end = [-1.;-1.]

      # initialize the dataframes to save the data
      data = DataFrame(iteration = zeros(Int16,length(grid)),
                  grid = grid,
                  Vint = 4π*Na.*(phi.^2 ./ grid.^2),
                  phi = zeros(Float64,length(grid)))

      energies = DataFrame(e1s = Float64[],
                        E1 = Float64[],
                        E2 = Float64[])

      @printf("\n\nStarting KS algorithm with Na = %0.4f, α = %0.2f\n", Na, α)

      # Start of the self-consistent iteration
      for t = 1:max_iter

            @printf("\nITERATION %d, Na = %0.4f\n", t, Na)

            Vint = 4π*Na.*(phi.^2 ./ ((grid).^2))# .+ 1e-11))
            Vgp = Vext .+ Vint

            # calculation of the eigenfunctions with the effective mean-field potential
            eigv_new, phi_new = Numerov(0, 1, grid, Vgp, bc_0=bc_0[1:2], bc_end=bc_end[1:2], Estep=Estep, verbose=verbose)

            # mixing of the new density with the old one
            phi_new =  α.*phi_new[:,1] .+ (1-α).*phi

            # next boundary conditions will be based on the current phi, mixed with the old
            # NOTE conditions at the end can be either these or a simple decaying exponential, they should both work
            bc_0 .=  α.*phi_new[1:2] .+ (1-α).*bc_0
            @show bc_end .=  α.*phi_new[end-1:end] .+ (1-α).*bc_end
            #bc_end .=  [-1.;-1.;-1.;-1.;-1.;-1.] # use default exponential in Numerov

            # Consistency check through the energy

            E1 = eigv_new[1]
            @show T_S(grid,phi)
            @show E_ext(grid,phi,Vext)
            @show E2 = T_S(grid,phi)

            verbose && @printf("E_eig = %f\tE_func = %f\tdiff = %f\n", E1, E2, E2-E1)

            # Check on the convergence by looking at how different is the new density
            @show delta = norm(phi_new .- phi)

            # save the found density function to be used and compared in the next iteration
            phi = phi_new

            # save partial results to data
            if (t%stride == 1) || stride == 1
                  verbose && println("Saving calculations to memory...\n")
                  data_tmp = DataFrame(iteration = t .* ones(Int16,length(grid)),
                              grid = grid,
                              Vint = 4π*Na.*(phi.^2 ./ grid.^2),
                              phi = phi)

                  append!(data,data_tmp)
            end
            push!(energies, [eigv_new; E1; E2])

            if delta < 1e-5
                  @printf("\nConvergence reached after %d steps with δ = %0.9f\n", t, delta)
                  break
            end
            if t == max_iter
                  @warn("\nExiting from the loop after $max_iter iterations, delta = $delta\n")
            end
      end

      return data, energies
end



rmax = 12
h = 2e-4
grid = Vector(h:h:rmax)
α = 0.1 # mixing coefficient of the densities
Vext = 0.5 .* (grid).^2

Threads.@threads for Na ∈ [0.01; 0.1; 1; 10; 100]
      @time data, energy = solve_GP(Na, α, grid, Vext; max_iter=280, stride=1, verbose=false)
      CSV.write("./Data/gpfunctions_$(α)_$Na.csv", data)
      CSV.write("./Data/gpenergy_$(α)_$Na.csv", energy)
end

#Juno.@profiler solve_GP(0.1, α, grid, Vext; max_iter=50, stride=2, verbose=false)
