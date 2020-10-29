using basicDFT, Plots, LinearAlgebra, Printf, DataFrames, CSV
include("functionals.jl")

# TODO
# better fix for rho_old
# fix delta
# check numerical instability in the flex of Vh
# fix energy calculation
# compare results with different initial condition

function solve_KS(N, rs, α, grid; max_iter=20, stride=2, verbose=false)

      Vext = V_ext.(grid,Rc(N,rs),rho_b(rs))
      Vext .= Vext #.+ abs(minimum(Vext))

      # set the initial trial electron density
      cos_single(x,l,c) = cos((x-c)*pi/l)^2 * (abs((x-c)*pi/l)<pi/2) + 1e-12
      rho = cos_single.(grid, 18, -1)
      rho = rho .* N ./ norm(rho,1)
      rho_old = rho

      # initial boundary conditions for the wavefunctions (s,p,d)
      # the bc for different values of l will become different, so here concatenate them
      bc_0 = [rho[1:2]./N; rho[1:2]./N; rho[1:2]./N].^(0.5)
      bc_end = -1 .* ones(Float64,3*2)

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

      # Start of the self-consistent Kohn-Sham method
      for t = 1:max_iter

            @printf("\nITERATION %d\n", t)

            data_step, energy_tmp = kohn_sham_step!(grid, Vext, rho, bc_0, bc_end, verbose=verbose)

            # Check on the convergence by looking at how different is the new density
            @show all(rho .== rho_old)
            @show delta = norm(rho .- rho_old)

            # save partial results to data
            replace!(data_step.iteration, -1 => Int16(t))
            t%stride==1 && append!(data,data_step)
            push!(energies, energy_tmp)

            # next boundary conditions will be based on the current eigenfunctions
            bc_0_new = [data_step.eigf_1s[1:2,1]; data_step.eigf_1p[1:2,1]; data_step.eigf_1d[1:2,1]]
            bc_end_new = [data_step.eigf_1s[end-1:end,1]; data_step.eigf_1p[end-1:end,1]; data_step.eigf_1d[end-1:end,1]]

            # mixing of the solution with the old one
            rho .=  α.*rho .+ (1-α).*rho_old
            @show bc_0 .=  α.*bc_0_new .+ (1-α).*bc_0
            @show bc_end .=  [-1.;-1.;-1.;-1.;-1.;-1.] #NOTE was needed to keep Numerov in check

            rho_old = rho # save the current density function for later mixing

            if (delta < 1e-8) && (t > 1)
                  @printf("\nConvergence reached after %d steps with δ = %f\n", t, delta)
                  #break
            end
      end
      return data, energies
end


# computes the mean-field potential Vks, solves the Shrodinger equation for the
# relevant quantum numbers, computes the new electronic density and returns everything
# as a dataframe
function kohn_sham_step!(grid::Vector, Vext::Vector, rho::Vector, bc_0::Vector, bc_end::Vector; Vks_cutoff=1e4, Estep=3e-3, verbose=false)

      #@show minimum(rho[findall(rho.>0.0)])
      # Hartree potential term, hotspot of the code
      Vh = V_h(grid, rho)
      # Kohn-Sham potential, function of rho
      Vks = Vext .+ Vh .+ V_xc(rho)
      # sharp cutoff on the potential
      Vks[Vks.>Vks_cutoff] .= Vks_cutoff

      # calculation of the eigenfunctions with the effective mean-field potential
      eigv_l0, eigf_l0 = Numerov(0, 2, grid, Vks, bc_0=bc_0[1:2], bc_end=bc_end[1:2], Estep=Estep, verbose=verbose)
      eigv_l1, eigf_l1 = Numerov(1, 1, grid, Vks, bc_0=bc_0[3:4], bc_end=bc_end[3:4], Estep=Estep, verbose=verbose)
      eigv_l2, eigf_l2 = Numerov(2, 1, grid, Vks, bc_0=bc_0[5:6], bc_end=bc_end[5:6], Estep=Estep, verbose=verbose)

      # compute the total electron density
      rho .= 2 .* eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1].^2 #.+ 1e-15
      rho .+= 2 .* eigf_l0[:,2].^2 + 10 .* eigf_l2[:,1].^2

      # save the computed functions (note that the vanilla, unmixed rho is saved here)
      data_tmp = DataFrame(iteration = -1 .* ones(Int16,length(grid)),
                        grid = grid,
                        Vh = Vh,
                        Vks = Vks,
                        rho = rho,
                        eigf_1s = eigf_l0[:,1],
                        eigf_2s = eigf_l0[:,2],
                        eigf_1p = eigf_l1[:,1],
                        eigf_1d = eigf_l2[:,1])

      # Consistency check through the energy
      E1 = E_ks(grid, rho, Vext, Vh)
      # sum of the eigenvalues - 1/2 hartree energy - exchange
      E2 = eigv_l0[1]*2 + eigv_l1[1]*6 - E_H(grid,rho,Vh) + E_XC(grid,rho)
      #verbose && @printf("E_ks = %f\tE_eig = %f\tdiff = %f\n", E1, E2, E2-E1)

      return data_tmp, [eigv_l0; eigv_l1[1]; eigv_l2[1]; E1; E2]
end


function V_ext(r::Float64, Rc::Float64, rho_b::Float64)
    if r > Rc
        return -4*pi*rho_b*Rc^3/(3*r)
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc)
    end
end


Rc(N,rs) = cbrt(N)*rs # radius of the positive jellium
rho_b(rs) = 3/(4π*rs^3) # density of charge inside the nucleus
N = 20
rs_Na = 3.93
rs_K = 4.86
rmax = 26
h = 5e-4
grid = Vector(h:h:rmax)
α = 0.2 # mixing coefficient of the densities

# Juno.@profiler
@time data, energy = solve_KS(N, rs_Na, α, grid, max_iter=20, stride=2, verbose=false)
CSV.write("./Data/ksfunctions_$N.csv", data)
CSV.write("./Data/ksenergy_$N.csv", energy)
