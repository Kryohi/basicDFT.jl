include("functionals.jl")
include("numerov.jl")


function solve_GP(Na, α, grid, Vext; max_iter=80, stride=1, verbose=false)

      # set the initial trial wavefunction
      cos_single(x,l,c) = cos((x-c)*π/l)^2 * (abs((x-c)*π/l) < π/2) + 1e-12
      phi = cos_single.(grid, 24, -1)
      phi = phi ./ norm(rho,1)

      # initial boundary condition for the wavefunction
      bc_0 = phi[1:2]
      bc_end = [-1.;-1.]

      # initialize the dataframes to save the data
      data = DataFrame(iteration = zeros(Int16,length(grid)),
                  grid = grid,
                  phi = zeros(Float64,length(grid)))

      energies = DataFrame(e1s = Float64[],
                        E1 = Float64[],
                        E2 = Float64[])

      @printf("\n\nStarting KS algorithm with Na = %d, α = %0.2f\n", N, α)

      # Start of the self-consistent iteration
      for t = 1:max_iter

            @printf("\nITERATION %d\n", t)

            Vgp = Vext .+ 4π*Na.*(phi.^2 ./ grid.^2)

            # calculation of the eigenfunctions with the effective mean-field potential
            eigv_new, phi_new = Numerov(0, 1, grid, Vgp, bc_0=bc_0[1:2], bc_end=bc_end[1:2], Estep=Estep, verbose=verbose)

            # mixing of the new density with the old one
            phi_new =  α.*phi_new .+ (1-α).*phi

            # next boundary conditions will be based on the current phi, mixed with the old
            # NOTE conditions at the end can be either these or a simple decaying exponential, they should both work
            bc_0 .=  α.*phi_new[1:2,1] .+ (1-α).*bc_0
            @show bc_end .=  α.*phi_new[end-1:end,1] .+ (1-α).*bc_end
            #bc_end .=  [-1.;-1.;-1.;-1.;-1.;-1.] # use default exponential in Numerov

            # Consistency check through the energy
            @show T_S(grid,phi)
            @show E_ext(grid,phi,Vext)
            @show E1 = T_S(grid,phi)

            E2 = eigv_l0[1]*2 + eigv_l1[1]*6 - E_H(grid,rho,Vh) - E_XC(grid,rho)

            verbose && @printf("E_ks = %f\tE_eig = %f\tdiff = %f\n", E1, E2, E2-E1)

            # Check on the convergence by looking at how different is the new density
            @show delta = norm(data_step.rho .- rho)

            # save the found density function to be used and compared in the next iteration
            phi = phi_new

            # save partial results to data
            replace!(data_step.iteration, -1 => Int16(t))
            t%stride==1 && append!(data,data_step)
            push!(energies, energy_step)

            if delta < 1e-6
                  @printf("\nConvergence reached after %d steps with δ = %0.9f\n", t, delta)
                  break
            end
            if t == max_iter
                  @warn("\nExiting from the loop after $max_iter iterations, delta = $delta\n")
            end
      end

      return data, energies
end



Na_ex = 10 .^ [-2:2]
rmax = 28
h = 5e-4
grid = Vector(h:h:rmax)
α = 0.2 # mixing coefficient of the densities
