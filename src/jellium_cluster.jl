using CSV
include("kohn_sham.jl")

Rc(N,rs) = cbrt(N)*rs # radius of the positive jellium
rho_b(rs) = 3/(4π*rs^3) # density of charge inside the nucleus

function V_ext(r::Float64, Rc::Float64, rho_b::Float64)
    if r > Rc
        return -4*pi*rho_b*Rc^3/(3*r)
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc)
    end
end

N = 20
rs_Na = 3.93
rs_K = 4.86
rmax = 28
h = 2.5e-4
grid = Vector(h:h:rmax)
α = 0.2 # mixing coefficient of the densities


Vext = V_ext.(grid, Rc(N,rs_Na), rho_b(rs_Na))
@time data, energy = solve_KS(N, α, grid, Vext, max_iter=60, stride=2)
CSV.write("./Data/ksfunctions_Na_$N.csv", data)
CSV.write("./Data/ksenergy_Na_$N.csv", energy)

Vext = V_ext.(grid, Rc(N,rs_K), rho_b(rs_K))
@time data, energy = solve_KS(N, α, grid, Vext, max_iter=60, stride=2)
CSV.write("./Data/ksfunctions_K_$N.csv", data)
CSV.write("./Data/ksenergy_K_$N.csv", energy)


# Juno.@profiler solve_KS(20, rs_K, α, grid, Vext, max_iter=10, stride=2)

#@code_native V_h(grid, last_rho)
#Juno.@profiler Vhtest(grid, rho)
#@btime Vhtest(grid, rho)
