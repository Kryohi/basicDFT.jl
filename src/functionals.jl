include("common_math.jl")

# TODO
# check everything
# more general handling of V_ext and local_energy
# always use simpson_integral()?

# Kohn-Sham potential
function V_ks(r::Float64, V_ext::Function, rho::Vector)

    # Exchange-correlation potential
    rho_r = rho[round(Int,h*r)]
    V_xc = (-0.25*cbrt(3/pi)*pow(rho_r,-2/3) - 0.14667*cbrt(4/(3*pi))*pow(rho_r,-2/3) / (0.78+3/(4*pi*rho_r))) * local_energy(rho_r) + local_energy(rho_r)

    return V_ext(r) + V_h(r, rho) + V_xc
end


# Hartree potential in the Kohn-Sham equation, also used to compute the Hartree energy
function V_h(r::Float64, rho::Vector)

    V_h = 0.
    x = round(Int,h*r)

    # Simpson integral
    for i = 2:2:x-1
        V_h += h * (rho[i-1]*(h*(i-1)) + 4*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3
    end

    for i = x:2:lenght(rho)-1
        V_h += h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / (3*r)
    end
    return 4*pi*V_h
end

function V_ext(r::Float64)

    if (r>Rc)
        return -4*pi*rho_b*Rc*Rc*Rc/(3*r)
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc)
    end
end


# Energy functional
function E_ks(rho::Vector)

    return T_S(rho) + E_ext(rho) + E_H(rho) + E_XC(rho)
end

# External energy
function E_ext(rho::Vector, V_ext::Function)

    rhoVext = rho .* V_ext.(0:h:h*length(rho)-h)

    return simpson_integral(rhoVext, length(rho), h);
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
