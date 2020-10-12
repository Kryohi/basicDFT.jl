using LinearAlgebra, Printf
include("common_math.jl")

const h2m = 0.5


# Performs the whole algorithm and finds the spectrum up to the n-th level
# Returns the found eingenvalues and eigenvectors

function Numerov(l, nmax, grid, V::Array; bc_0=[-1.,-1.], bc_end=[-1.,-1.], Estep=-1., verbose=false)

    if length(grid) != length(V)
        throw(ArgumentError(grid, "length of grid and V must match"))
    end
    if (length(bc_0) != 2)
        throw(ArgumentError(bc_0, "boundary conditions must contain two points"))
    end
    if any(x->isnan(x), V)
        throw(error("Potential V contains NaNs"))
    end

    xmin = length(bc_0)+1
    xmax = length(grid)
    h = grid[2]-grid[1]

    # centrifugal term, dependent on l but not on E
    centrifugal = zeros(Float64, xmax)
    # kinetic term
    k2 = zeros(Float64, xmax)
    # calculated wavefunction, respectively before xc and after xc (f stands for forward)
    yf = zeros(Float64, xmax)
    yb = zeros(Float64, xmax)
    # eigenfunctions found
    eigf = zeros(Float64, xmax, nmax)
    # eigenvalues found
    eigv = zeros(Float64, nmax)

    # precalculation of the known terms (independent of E)
    for x = 1:xmax
        centrifugal[x] = l*(l+1)/(x*h*x*h+1e-14)
    end

    # boundary conditions, if provided (otherwise exponential decay will be used)
    bc_0_exp = (bc_0[1] == -1.)
    bc_end_exp = (bc_end[1] == -1.)
    if !bc_0_exp
        yf[1] = bc_0[1]
        yf[2] = bc_0[2]
    end
    if !bc_end_exp
        yb[xmax-1] = bc_end[1] # check if the order is correct
        yb[xmax] = bc_end[2]
    end
    # else they have to be computed at each step (they depend on E)

    # if Estep is not provided as argument, we use the small difference of V near the minimum
    Vmin, Vmin_idx = findmin(V)
    println(Vmin)
    println(Vmin_idx)
    (Estep == -1.0) && (Estep = abs(Vmin-V[Vmin_idx+1])*100)
    # starting (inferior) energy for the Numerov algorithm
    E = Vmin+Estep+1e-9
    verbose && @printf("Starting at E = %f", E)
    # number of eigenvalues found
    nfound = 1
    # temporary variable with the previous value of the log derivative difference of yf and yb
    prevdelta = 0.
    # needed because when function flips sign the delta also changes sign producing false positives
    prev_yc = 0.

    # Iterative process to find the energy eigenvalues and the corresponding eigenvectors
    while (nfound <= nmax)

        # we search for the point of intersection of V with the current E
        _, xc = findmin(abs.(reverse(E .- V)))
        xc = xmax - xc
        (xc > length(V)-3) && (throw(error("xc at the end of the domain")))
        #xc, _ = secant(V .- E, (xmax*9)÷10, xmax-1, 10, 10^3)
        #xc, _ = secant(V .- E, 1, xmax-1, 10, 10^3)

        # we run forward and backward Numerov, store the results in yf and yb
        # and compute the difference in the derivative of the logarithms at xc
        delta = findDelta!(E, V, centrifugal, k2, h, xmin, xmax, bc_0_exp, bc_end_exp, verbose, yf, yb)

        # We are searching for the zeros of delta, with increasing values of E, so
        # if there is a change in sign, start the finer search of the 0 of delta(E)
        # with the secant method, between E and E-Estep
        if (delta*prevdelta < 0 && yf[xc]*prev_yc > 0)
            verbose && @printf("\n[Numerov] Found a point of inversion at %0.9f - %0.9f\n", E-Estep, E)

            eigv[nfound], _ = secant(e ->
                findDelta!(e, V, centrifugal, k2, h, xmin, xmax, bc_0_exp, bc_end_exp, verbose, yf, yb),
                E-Estep, E, 1e-5, 10^4)


            verbose && @printf("\nE%d = %.9f\n\n", nfound, eigv[nfound])
            _, xc = findmin(abs.(eigv[nfound] .- V))
            (xc > length(V)-3) && (throw(error("xc at the end of the domain")))

            # put together yf and yb to form the eigenfunction
            for x = 1:xc
                eigf[x,nfound] = yf[x]
            end
            for x = xc+1:xmax
                eigf[x,nfound] = yb[x] * yf[xc]/yb[xc]  #impose continuity
            end

            # normalize the eigenfunction
            norm2 = norm(eigf[:,nfound])
            verbose && println("norm=",norm2)
            eigf[:,nfound] = eigf[:,nfound] ./ norm2

            # update the number of solutions found
            nfound += 1
        end

        prevdelta = delta
        prev_yc = yf[xc]
        E += Estep

        if (E > maximum(V))
           error("[Numerov] Could not find enough solutions with the parameters provided\n")
        end
    end

    return eigv, eigf
end


function findDelta!(E, V::Vector, centrifugal, k2, h, xmin, xmax, bc_0_exp, bc_end_exp, verbose, yf, yb)

    # Boundary conditions at rmin
    if bc_0_exp
        yf[1] = exp(-sqrt(abs(E)/h2m)*xmax*h)
        yf[2] = exp(-sqrt(abs(E)/h2m)*(xmax-1)*h)
    end
    # Boundary conditions at rmax
    if bc_end_exp
        yb[xmax] = exp(-sqrt(abs(E)/h2m)*xmax*h)
        yb[xmax-1] = exp(-sqrt(abs(E)/h2m)*(xmax-1)*h)
    end

    # we search for the last point of intersection of V with the current E
    _, xc = findmin(abs.(reverse(E .- V)))
    xc = xmax - xc
    (xc > length(V)-3) && (throw(error("xc at the end of the domain")))
    verbose && println("xc = ", xc)
    #xc, _ = secant(V .- E, (xmax*9)÷10, xmax-1, 10, 10^3)

    # calculation of the
    k2 .= (E .- V) ./ h2m .- centrifugal

    numerov_forward!(h, xc, xmin, k2, yf)
    numerov_backward!(h, xc, xmax, k2, yb)

    delta = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc]

    verbose && @printf("yf[xc-1] = %0.9f,\tyf[xc] = %0.9f,\tyf[xc+1]=%0.9f\n", yf[xc-1], yf[xc], yf[xc+1])
    verbose && @printf("yb[xc-1] = %0.9f,\tyb[xc] = %0.9f,\tyb[xc+1]=%0.9f\n", yb[xc-1], yb[xc], yb[xc+1])
    verbose && @printf("derforward = %0.12f,  derback = %0.12f\n", der5(yf,xc,h)/yf[xc],der5(yb,xc,h)/yb[xc])
    verbose && @printf("E = %f\tDelta rough = %f\n", E, delta)

    return delta
end


# finds the
# xmin is used in order to not overwrite yf at the boundary provided to Numerov
function numerov_forward!(h::Float64, xc::Int64, xmin::Int64, k2, yf)

    hh = h*h
    c0 = hh*k2[xmin]/12
    c_1 = hh*k2[xmin-1]/12
    c_2 = hh*k2[xmin-2]/12

    @inbounds for x = xmin:xc+3
        yf[x] = (yf[x-1]*(2-10*c_1) - yf[x-2]*(1+c_2)) / (1+c0)
        c_2 = c_1
        c_1 = c0
        c0 = hh*k2[x+1]/12
    end
end

function numerov_backward!(h::Float64, xc::Int64, xmax::Int64, k2, yb)

    hh = h*h
    c0 = hh*k2[xmax-2]/12
    c1 = hh*k2[xmax-1]/12
    c2 = hh*k2[xmax]/12

    @inbounds for x = xmax-2:-1:xc-3
        yb[x] = (yb[x+1]*(2-10*c1) - yb[x+2]*(1+c2)) / (1+c0)
        #(x < xc+200) && println("yb = ", yb[x])
        c2 = c1
        c1 = c0
        c0 = hh*k2[x-1]/12
    end
end


# TODO:
# understand if it is ever useful to have longer than 2 boundary conditions, simplify code
# better explain the need for yc
# check boundary conditions for 1D HO
