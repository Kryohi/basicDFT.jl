include("common_math.jl")

# TODO
# check everything
# more general handling of V_ext and local_energy
# always use simpson_integral()?
# 0.78 or 7.8?


# Local energy used in the LDA
@inline function local_energy(rho::Float64)
    return -0.75*cbrt(3*rho/pi) - 0.44/(7.8 + 3/(4*pi*rho)) #TODO check if ∛
end

# Exchange-Correlation potential
function V_xc(rho::Vector)
    # derivative of xc energy times rho
    dExc = -0.25*cbrt(3/π).*cbrt.(rho) - 0.44*2^(2/3) ./ (3 .*cbrt.(3π.*rho).*(7.8 .+ 2^(2/3)./cbrt.(3π.*rho)).^2)

    return local_energy.(rho) .+ dExc
end


# Hartree potential in the Kohn-Sham equation, also used to compute the Hartree energy

@inbounds function V_h(grid::Vector, rho::Vector)
    Vh = zeros(Float64, length(rho))
    h = grid[2]-grid[1]

    Vh[end] = 4pi*simpson_integral(grid.^2 .* rho, h)/grid[1]
    Vh[1] = 4pi*simpson_integral(rho.*grid, h)

    # this is the most expensive function to compute, so we use multithreading
    # (which unexpectedly seems to work without much user input, with a 40% faster program using 2 cores)
    Threads.@threads for x = 2:length(grid)-1
        for i = 2:2:x-1
            Vh[x] += h * (rho[i-1]*h*(i-1)*h*(i-1) + 4.0*rho[i]*h*i*h*i + rho[i+1]*h*(i+1)*h*(i+1)) / (3.0* grid[x])
        end
        for i = x:2:length(rho)-1
            Vh[x] += h * (rho[i-1]*h*(i-1) + 4.0* rho[i]*h*i + rho[i+1]*h*(i+1)) / 3.
        end
        Vh[x] = 4pi*Vh[x]  # from spherical integration
    end

    return Vh
end

# # Hartree potential in the Kohn-Sham equation, also used to compute the Hartree energy - alternative version, may be faster
# @inbounds function V_h(grid::Vector, rho::Vector)
#     Vh = zeros(Float64, length(rho))
#     Vh1 = zeros(Float64, length(rho))
#     Vh2 = zeros(Float64, length(rho))
#     h = grid[2]-grid[1]
#
#     for i = 2:2:length(grid)-1
#         Vh1[i] = h * (rho[i-1]*(h*(i-1)) + 4*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3
#     end
#     for i = 2:2:length(grid)-1
#         Vh2[i] = h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / 3
#     end
#     Vh1[1] = h*rho[1]*grid[1]
#     Vh1[end] = Vh1[end] + h*rho[end]*grid[end]
#     Vh2[1] = h*rho[1]*grid[1]^2
#     Vh2[end] = Vh1[end] + h*rho[end]*grid[end]^2
#
#     Vh1 .= (4*pi) .* Vh1
#     Vh2 .= (4*pi) .* Vh2
#
#     Vh[1] = sum(Vh2)/grid[1]
#     for i = 2:length(grid)
#         Vh[i] = sum(Vh1[1:i-1]) + sum(Vh2[i:length(grid)])/grid[i]
#     end
#     #println("Vh[1] = ", Vh[1])
#     return Vh
# end



# Energy functional
function E_ks(grid::Vector, rho::Vector, Vext::Vector, Vh::Vector)

    return T_S(grid,rho) + E_ext(grid,rho, Vext) + E_H(grid,rho,Vh) + E_XC(grid,rho)
end

# External energy
function E_ext(grid::Vector, rho::Vector, V_ext::Vector)
    h = grid[2]-grid[1]
    return simpson_integral(rho .* V_ext, length(rho), h)
end

# Hartree energy TODO check it bc i'm really not sure why i'm programming at BUC 10 minutes before closure
function E_H(grid::Vector, rho::Vector)

    h = grid[2]-grid[1]
    E_h = simpson_integral(rho .* V_h(grid, rho), h) / 2

    return E_h  #check this
end
# same function but reusing V_H precalculated
function E_H(grid::Vector, rho::Vector, Vh::Vector)

    h = grid[2]-grid[1]
    E_h = simpson_integral(rho .* Vh, h) / 2

    return E_h  #check this
end


# Exchange-correlation energy
function E_XC(grid::Vector, rho::Vector)

    h = grid[2]-grid[1]
    E_xc = simpson_integral(rho .* local_energy.(rho), h)

    return E_xc
end

# Kohn–Sham kinetic energy (NOTE should take the ks orbitals as input?)
function T_S(grid::Vector, phi::Vector)
    gridlength = length(phi)
    h = grid[2]-grid[1]
    integrand = der5(phi,h).^2
    h2m = .5
    return h2m * simpson_integral(integrand, gridlength, h);
end
