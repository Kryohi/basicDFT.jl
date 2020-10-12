include("common_math.jl")

# TODO
# check everything
# more general handling of V_ext and local_energy
# always use simpson_integral()?

# Kohn-Sham potential

function V_ks(grid::Vector, V_ext::Vector, rho::Vector)
    h = grid[2]-grid[1]

    # Exchange-correlation potential
    # replaces all zeros in rho
    rho[findall(rho.==0.0)] .= minimum(rho[findall(rho.>0.0)])
    rho23 = rho .^ (-2/3)
    println("rho23[1] = ", rho23[1])
    V_xc = (-0.25*cbrt(3/pi) .* rho23 .- 0.14667 .* cbrt.(4/(3*pi)) .* rho23 ./ (0.78 .+ 3 ./ (4*pi.*(rho.+1e-12)))) .* local_energy.(rho) .+ local_energy.(rho)

    Vks = V_ext .+ V_h(grid, rho) .+ V_xc
    println("Vxc[1] = ", V_xc[1])

    return Vks
end


# Hartree potential in the Kohn-Sham equation, also used to compute the Hartree energy
#
# function V_h(grid::Array, rho::Vector)
#     Vh = zeros(Float64, length(rho))
#     h = grid[2]-grid[1]
#
#     for i = 2:2:length(grid)-1
#         Vh[end] += h * (rho[i-1]*(h*(i-1)) + 4*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3
#     end
#     for i = 2:2:length(grid)-1
#         Vh[1] += h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / (3*grid[2])
#     end
#     Vh[1] = 4*pi*Vh[1]
#
#     for x = 2:length(grid)-1
#         # Simpson integral
#         for i = 2:2:x-1
#             Vh[x] += h * (rho[i-1]*(h*(i-1)) + 4*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3
#         end
#         for i = x:2:length(rho)-1
#             Vh[x] += h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / (3*grid[x+1]) # what do in 0
#         end
#         Vh[x] = 4*pi*Vh[x]
#     end
#
#     return Vh
# end


function V_h(grid::Array, rho::Vector)
    Vh = zeros(Float64, length(rho))
    Vh1 = zeros(Float64, length(rho))
    Vh2 = zeros(Float64, length(rho))
    h = grid[2]-grid[1]

    for i = 2:2:length(grid)-1
        Vh1[i] = h * (rho[i-1]*(h*(i-1)) + 4*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3
    end
    for i = 2:2:length(grid)-1
        Vh2[i] = h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / 3
    end
    Vh1[1] = h*rho[1]*grid[1]
    Vh1[end] = Vh1[end] + h*rho[end]*grid[end]
    Vh2[1] = h*rho[1]*grid[1]^2
    Vh2[end] = Vh1[end] + h*rho[end]*grid[end]^2

    Vh1 .= (4*pi) .* Vh1
    Vh2 .= (4*pi) .* Vh2

    Vh[1] = sum(Vh2)/grid[2]
    for i = 2:length(grid)
        Vh[i] = sum(Vh1[1:i-1]) + sum(Vh2[i:length(grid)])/grid[i]
    end
    println("Vh[1] = ", Vh[1])
    return Vh
end



# Energy functional
function E_ks(rho::Vector, V_ext::Vector)

    return T_S(rho) + E_ext(rho, V_ext) + E_H(rho) + E_XC(rho)
end

# External energy
function E_ext(rho::Vector, V_ext::Vector)
    return simpson_integral(rho .* V_ext, length(rho), h)
end

# Hartree energy TODO check it bc i'm really not sure why i'm programming at BUC 10 minutes before closure
function E_H(rho::Vector)

    E_h = 0.
    integrand = rho .* V_h.(0:h:h*length(rho)-h)

    for i = 2:2:length(rho)-1
        E_h += h*(rho[i-1]*integrand[i-1] + 4*rho[i]*integrand[i] + rho[i+1]*integrand[i+1])/3
    end

    return E_h/2  #check this
end

# Exchange-correlation energy
function E_XC(rho::Vector)

    E_xc = 0.
    # simpson integration
    for i = 2:2:length(rho)-1
        E_xc += h*(rho[i-1]*local_energy(rho[i-1]) + 4*rho[i]*local_energy(rho[i]) + rho[i+1]*local_energy(rho[i+1]))/3
    end

    return E_xc;
end


# Kohn–Sham kinetic energy (should take the ks orbitals as input)
function T_S(phi::Vector)

    #excludes extrema used as boundary conditions
    integrand = zeros(Float64, gridlength-4)

    for i = 3:gridlength-2
        integrand[i-2] = der5(phi,i,h) * der5(phi,i,h)
    end

    return h2m * simpson_integral(integrand, gridlength-4, h);
end


# Local energy used in the LDA
@inline function local_energy(rho::Float64)
    return -0.75*cbrt(3*rho/pi) - 0.44/(0.78 + 3/(4*pi*rho))
end
