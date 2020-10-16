using basicDFT, Plots, LinearAlgebra, Printf, DataFrames, CSV
include("functionals.jl")

# TODO
# capire perche bcend piccolo rompe numerov?
# numerov from 0 to Vks cutoff?
# fix multiple numerov solutions
# add 2s orbital
# understand and implement N=20 case

N = 8
rs_Na = 3.93
Rc(rs) = cbrt(N)*rs # radius of the positive jellium
rho_b(rs) = 3/(4*pi*rs^3) # density of charge inside the nucleus
rmax = 20
h = 5e-4
grid = Vector(h:h:rmax)
α = 0.2 # mixing coefficient of the densities


function solve_KS(N, rs, α, grid; max_iter=20, verbose=false)

      Vext = V_ext.(grid,Rc(rs),rho_b(rs))
      Vext .= Vext #.+ abs(minimum(Vext))

      # set the initial trial electron density and boundary conditions
      cos_single(x,l,c) = cos((x-c)*pi/l)^2 * (abs((x-c)*pi/l)<pi/2) + 1e-9
      rho = cos_single.(grid, 12, Rc(rs))
      rho = rho .* 8 ./ norm(rho,1)

      # initial boundary conditions for the wavefunctions
      # the bc for different values of l will become different, so here concatenate them
      bc_0 = [rho[1:2]./6; rho[1:2]./6]
      bc_end = -1 .* ones(Float64,4)

      # initialize the dataframe to save the data
      data = DataFrame(iteration = zeros(Int16,length(grid)),
                  grid = grid,
                  Vks = zeros(Float64,length(grid)),
                  rho = rho,
                  eigf_1s = zeros(Float64,length(grid)),
                  eigf_1p = zeros(Float64,length(grid)));

      # Start of the self-consistent Kohn-Sham method
      for t = 1:max_iter

            @printf("\nITERATION %d\n", t)
            rho_old = rho # save the current density function for later mixing

            data_step = kohn_sham_step!(grid, Vext, rho, bc_0, bc_end, verbose=verbose)

            # Check on the convergence by looking at how different is the new density
            @show delta = norm(rho .- rho_old)

            # save partial results to data
            replace!(data_step.iteration, -1 => Int16(t))
            append!(data,data_step)

            # next boundary conditions will be based on the current eigenfunctions
            bc_0_new = [data_step.eigf_1s[1:2,1]; data_step.eigf_1p[1:2,1]]
            bc_end_new = [data_step.eigf_1s[end-1:end,1]; data_step.eigf_1p[end-1:end,1]]

            # mixing of the solution with the old one
            rho .=  α.*rho .+ (1-α).*rho_old
            @show bc_0 .=  α.*bc_0_new .+ (1-α).*bc_0
            @show bc_end .=  [-1.;-1.;-1.;-1.] #NOTE was needed to keep the Numerov in check

            if delta < 1e-6
                  @printf("\nConvergence reached after %d steps with δ = %f\n", t, delta)
                  #break
            end
      end
      return data
end


# computes the mean-field potential Vks, solves the Shrodinger equation for the
# relevant quantum numbers, computes the new electronic density and returns everything
# as a dataframe
function kohn_sham_step!(grid::Vector, Vext::Vector, rho::Vector, bc_0::Vector, bc_end::Vector; Vks_cutoff=1e4, Estep=3e-3, verbose=false)

      # total Kohn-Sham potential, function of rho
      Vks = V_ks(grid, Vext, rho)
      # sharp cutoff on the potential
      Vks[Vks.>Vks_cutoff] .= Vks_cutoff

      # calculation of the eigenfunctions with the effective mean-field potential
      eigv_l0, eigf_l0 = Numerov(0, 2, grid, Vks, bc_0=bc_0[1:2], bc_end=bc_end[1:2], Estep=Estep, verbose=verbose)
      eigv_l1, eigf_l1 = Numerov(1, 2, grid, Vks, bc_0=bc_0[3:4], bc_end=bc_end[3:4], Estep=Estep, verbose=verbose)

      # compute the total electron density
      rho .= 2. *eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1].^2 .+ 1e-9#./ 3

      # save the computed functions (note that the vanilla, unmixed rho is saved here)
      data_tmp = DataFrame(iteration = -1 .* ones(Int16,length(grid)),
                        grid = grid,
                        Vks = Vks,
                        rho = rho,
                        eigf_1s = eigf_l0[:,1],
                        eigf_1p = eigf_l1[:,1])

      # Consistency check through the energy
      #E1 = E_ks(grid, rho, Vext)
      # sum of the eigenvalues - 1/2 hartree energy - exchange
      #E2 = eigv_l0[1]*2 + eigv_l1[1]*6 - E_H(grid,rho) + E_XC(grid,rho)
      #verbose && @printf("E_ks = %f\tE_eig = %f\tdiff = %f\n", E1, E2, E2-E1)

      return data_tmp
end


function V_ext(r::Float64, Rc::Float64, rho_b::Float64)
    if r > Rc
        return -4*pi*rho_b*Rc^3/(3*r)
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc)
    end
end

# Juno.@profiler
@time data = solve_KS(N, rs_Na, α, grid, max_iter=24)
CSV.write("./Data/ksfunctions.csv", data)
