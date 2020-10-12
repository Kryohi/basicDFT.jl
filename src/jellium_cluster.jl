using basicDFT, Plots, Printf
include("functionals.jl")

# TODO
# understand boundary conditions, parity, if to shift the potential, centrifugal force
# understand and implement N=20 case


function V_ext(r::Float64, Rc::Float64, rho_b::Float64)
    if (r>Rc)
        return -4*pi*rho_b*Rc*Rc*Rc/(3*r)
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc)
    end
end


N = 8
rs_Na = 3.93
Rc(rs) = cbrt(N)*rs
rho_b(rs) = 3/(4*pi*rs*rs*rs) # density of charge inside the nucleus

rmax = 20
h = 5e-4
grid = h:h:rmax
rho = zeros(Float64, length(grid))
Vks = zeros(Float64, length(grid)) # total Kohn-Sham potential, function of rho

α = 0.1 # mixing coefficient of the densities
Vext = V_ext.(grid,Rc(rs_Na),rho_b(rs_Na))
Vext .= Vext .+ abs(minimum(Vext))


# starting guess for the density is calculated from a system in V_ext with no e-e interactions

_, eigf_l0 = Numerov(0, 2, grid, Vext, bc_0=[0.,1.], Estep=1e-2, verbose=false)
_, eigf_l1 = Numerov(1, 2, grid, Vext, bc_0=[0.,h], Estep=1e-2, verbose=false)

# calculate the initial electron density
rho .= 2. *eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1] ./ 3

# save the current density function for later mixing
rho_old = rho

# Normalization


# Start of the self-consistent Kohn-Sham method
for t = 1:100

      @printf("\nITERATION %d\n", t)

      Vks = V_ks(collect(grid), Vext, rho)
      println(Vks[1])
      plot(grid, Vks)
      eigv_l0, eigf_l0 = Numerov(0, 2, grid, Vks, bc_0=[0.,1.], Estep=1e-2, verbose=true)
      eigv_l1, eigf_l1 = Numerov(1, 2, grid, Vks, bc_0=[0.,h], Estep=1e-2, verbose=true)

      rho .= 2. *eigf_l0[:,1].^2 + 6 .* eigf_l1[:,1] ./ 3

      # mixing of the solution with the old density (va messa prima o dopo ks?)
      rho .=  α .* rho + (1-α) .* rho_old
      #fprintf(ksdensity, "%f\n", RHO[i]);
      #rho_old .= rho

      # Consistency check through the energy
      #E1 = E_ks(rho, Vext)
      # sum of the eigenvalues - 1/2 hartree energy - exchange
      #E2 = eigv_l0[1]*2 + eigv_l1[1]*6 - E_H(rho) + E_XC(rho)
      #@printf("E_ks = %f\tE_ = %f\tdiff = %f\n", E1, E2, E2-E1)

      # Check on the convergence by looking at how different is the new density
      delta = maximum(abs.(rho .- rho_old))
      #fprintf(ksdensity, "%f\n", delta)

      if delta < 1e-4
            @printf("\nConvergence reached after %d steps.", t)
            break
      end
end
