

# Kohn-Sham potential
function V_ks(double r, const double *rho)

    # Exchange-correlation potential
    double V_xc = (-0.25*pow(3/pi,1./3)*pow(rho[(int)(h*r)],-2./3) - 0.14667*pow(4/(3*pi),1./3)*pow(rho[(int)(h*r)],-2./3)/(0.78+3/(4*pi*rho[(int)(h*r)])))*local_energy(rho[(int)(h*r)]) + local_energy(rho[(int)(h*r)]);

    return V_ext(r) + V_h(r, rho) + V_xc;
end


function V_ext(double r)

    if (r>Rc)
        return -4*pi*rho_b*Rc*Rc*Rc/(3*r);
    else
        return 2*pi*rho_b*(r*r/3-Rc*Rc);
end

# Hartee potential in the Kohn-Sham equation, also used to compute the Hartree energy
function V_h(double r, const double *rho)

    double V_h = 0.;
    for (int i=1; i < (int)(h*r)-1; i+=2)
        V_h += h * (rho[i-1]*(h*(i-1)) + 4.*rho[i]*(h*i) + rho[i+1]*(h*(i+1))) / 3.;

    for (int i=(int)(h*r); i < gridlength-1; i+=2)
        V_h += h * (rho[i-1]*(h*(i-1))*(h*(i-1)) + 4.*rho[i]*(h*i)*(h*i) + rho[i+1]*(h*(i+1))*(h*(i+1))) / (3.*r);

    return 4*pi*V_h;
end



# Energy functional
function E_ks(double *rho)

    return T_S(rho) + E_ext(rho) + E_H(rho) + E_XC(rho);
end

# External energy
function E_ext(double *rho)

    double rhoVext[gridlength];
    for (int i=0; i < gridlength; i++)
        rhoVext[i] = rho[i]*V_ext(i*h)
    end

    return simpson_integral(rhoVext, gridlength, h);
end

# Hartree energy TODO check it bc i'm really not sure why i'm programming at BUC 10 minutes before closure
function E_H(double *rho)

    double E_h = 0;
    double integrand[gridlength];

    for (int i=0; i < gridlength; i++)
        integrand[i] = V_h(i*h, rho)*rho[i];

    for (int i=1; i < gridlength-1; i+=2)
        E_h += h*(rho[i-1]*integrand[i-1] + 4.*rho[i]*integrand[i] + rho[i+1]*integrand[i+1])/3.;

    return E_h/2; #check this
end

# Exchange-correlation energy
function E_XC(double *rho)

    double E_xc = 0;
    for (int i=1; i < gridlength-1; i+=2)
        E_xc += h*(rho[i-1]*local_energy(rho[i-1]) + 4.*rho[i]*local_energy(rho[i]) + rho[i+1]*local_energy(rho[i+1]))/3.;

    return E_xc;
end


# Kohn–Sham kinetic energy (should take the ks orbitals as input)
function T_S(double *phi)

    integrand[gridlength-4]; #excludes extrema used as boundary conditions

    for(int i=2; i<gridlength-2; i++)
        integrand[i-2] = der5(phi,i,h) * der5(phi,i,h);

    return h2m * simpson_integral(integrand, gridlength-4, h);
end


# Local energy used in the LDA
@inline function local_energy(rho)
    return -0.75*pow(3*rho/pi, 1/3) - 0.44/(0.78 + 3/(4*pi*rho));
end
