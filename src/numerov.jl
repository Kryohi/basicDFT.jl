using LinearAlgebra, Printf
include("common_math.jl")

const h2m = 0.5


# Performs the whole algorithm and finds the spectrum up to the n-th level
# Returns the found eingenvalues and eigenvectors

function Numerov(l::Int, nmax::Int, grid, V::Vector; bc_0=[-1.,-1.], bc_end=[-1.,-1.], Estep=-1., tol=5e-6, maxiter=10^5, strict_xc=false, verbose=false)

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
    # calculated wavefunction, respectively before xc and after xc (f stands for forward)
    yf = zeros(Float64, xmax)
    yb = zeros(Float64, xmax)
    # eigenfunctions found
    tmp_eigf = zeros(Float64, xmax)
    eigf = zeros(Float64, xmax, nmax)
    # eigenvalues found
    eigv = zeros(Float64, nmax)

    # precalculation of the known terms (independent of E)
    for x = 1:xmax
        centrifugal[x] = l*(l+1)/(x*h*x*h+1e-14)
    end

    # boundary conditions, if provided (otherwise exponential decay will be used)
    bc_0_exp = (bc_0[1] ≈ -1.)
    bc_end_exp = (bc_end[1] ≈ -1.)
    if !bc_0_exp
        yf[1] = bc_0[1]
        yf[2] = bc_0[2]
    end
    if !bc_end_exp
        yb[xmax-1] = bc_end[1]
        yb[xmax] = bc_end[2]
    end
    # else they have to be computed at each step (they depend on E)

    # if Estep is not provided as argument, we use the small difference of V near the minimum
    Vmin, Vmin_idx = findmin(V)

    (Estep == -1.0) && (Estep = abs(Vmin-V[Vmin_idx+1])*100)
    # starting (inferior) energy for the Numerov algorithm
    E = Vmin + Estep/10 # + Estep
    verbose && @printf("Starting at E = %f", E)
    # number of eigenvalues found
    nfound = 1
    # temporary variable with the previous value of the log derivative difference of yf and yb
    prevdelta = 0.
    # flag put to true if the change of sign of delta comes with high values,
    # and a finer search with reduced Estep is performed
    # if Estep was not reduced, the secant method might not converge
    high_derivative_workaround = false
    # flag put to true if the number of nodes of the solution can't reach the correct
    # number after some attempts, then the solution with the lowest nodes is chosen
    max_attempts = 6
    n_attempts = 0
    nodes_backup = ones
    tmp_eigf_backup = zeros(Float64, xmax)
    eigv_backup = 0.0


    # Iterative process to find the energy eigenvalues and the corresponding eigenvectors
    while (nfound <= nmax)

        # we search for the point of intersection of V with the current E
        xc = findEnergyIntersection(E,V,strict_xc)

        # we run forward and backward Numerov, store the results in yf and yb
        # and compute the difference in the derivative of the logarithms at xc
        delta = findDelta!(E, V, centrifugal, h, xmin, xmax, bc_0_exp, bc_end_exp, strict_xc, verbose, yf, yb)

        (delta*prevdelta < 0) && (@show delta, prevdelta, E)

        # if the change of sign of delta comes with high values, a finer search
        # in a smaller interval is performed
        # if Estep was not reduced, the secant method might not converge
        if ((delta*prevdelta < 0) && abs(delta-prevdelta) > 100.0)
            high_derivative_workaround = true
            Estep = Estep/10
            for i=10:-1:1
                delta = findDelta!(E-Estep*i, V, centrifugal, h, xmin, xmax, bc_0_exp, bc_end_exp, strict_xc, verbose, yf, yb)

                if delta*prevdelta < 0
                    @info "finer search due to high derivative finished with \n\tdelta=$delta at E=$(E-Estep*i), prevDelta=$prevdelta at E=$(E-Estep*(i+1))"
                    E = E-Estep*i
                    break
                else
                    prevdelta = delta
                end
            end
        end

        # We are searching for the zeros of delta, with increasing values of E, so
        # if there is a change in sign, start the finer search of the 0 of delta(E)
        # with the secant method, between E and E-Estep
        if delta*prevdelta < 0 #&& abs(delta-prevdelta) < 1000.0
            verbose && @printf("\n[Numerov] Found a point of inversion at %0.9f - %0.9f\n", E-Estep, E)

            eigv[nfound], _ = secant(e ->
                findDelta!(e, V, centrifugal, h, xmin, xmax, bc_0_exp, bc_end_exp, strict_xc, verbose, yf, yb),
                E-Estep, E, tol, maxiter) #tol+1e-1*(nfound+1)÷2

            verbose && @printf("\nE%d = %.9f\n\n", nfound, eigv[nfound])

            # sometimes the secant gives NaNs due to a bad energy interval,
            # but the search for solutions can continue
            if isnan(eigv[nfound])
                # if the secant gives NaN but we have a previous solution with the wrong
                # number of nodes, we accept that
                if n_attempts>max_attempts-1
                    nn = count_nodes(tmp_eigf_backup)
                    @warn "accepting solution with $nn nodes at E$nfound = $eigv_backup"
                    # update chosen solution
                    yf = tmp_eigf_backup
                    yb = tmp_eigf_backup
                    eigv[nfound] = eigv_backup
                    # attempts counter increased in order to accept and save the solution later
                    n_attempts = 1000
                else
                    @warn "last eigenfunction search gave NaNs, skipping this energy interval \n\t high_derivative_workaround=$high_derivative_workaround"
                    @goto skip_interval
                end
            end

            xc = findEnergyIntersection(eigv[nfound],V,strict_xc)

            # put together yf and yb to form the eigenfunction
            for x = 1:xc
                tmp_eigf[x] = yf[x]
            end
            for x = xc+1:xmax
                tmp_eigf[x] = yb[x] * yf[xc]/yb[xc]  #impose continuity
            end

            # normalize the eigenfunction
            norm2 = sqrt(simpson_integral(tmp_eigf.^2, h))
            tmp_eigf = tmp_eigf ./ norm2

            # we count the nodes of the function as an additional test
            n_nodes = count_nodes(tmp_eigf)

            if n_nodes != nfound-1 || n_attempts > max_attempts
                @info "number of nodes of E$nfound is $n_nodes, skipping solution at $n_attempts attempts, "

                if n_attempts == 0 # should also probably check if solution is identical
                    eigv_backup = eigv[nfound]
                    tmp_eigf_backup = tmp_eigf
                end
                n_attempts += 1
                #nodes_backup[n_attempts] = n_nodes
                # we accept and save the first solution even if it might be wrong
                if n_attempts > max_attempts+1
                    if high_derivative_workaround
                        Estep = Estep*10
                        high_derivative_workaround = false
                    end
                    nn = count_nodes(tmp_eigf_backup)
                    @warn "accepting solution with $nn nodes at E$nfound = $eigv_backup, after $n_attempts tries"
                    # update chosen solution
                    eigf[:,nfound] = tmp_eigf_backup
                    eigv[nfound] = eigv_backup
                    # reset attempts counter
                    n_attempts = 0
                    # update the number of solutions found
                    nfound += 1
                end

            elseif nfound>1 && (abs(eigv[nfound]-eigv[nfound-1])<1e-5)

                @info "skipping equivalent solution at $(eigv[nfound]), with $n_nodes nodes"
                if n_nodes == nfound-1
                    @warn "solution discarded, but it has the correct number of nodes"
                    @show tmp_eigf[xc-5:xc+5]
                end
            else
                # once in a while, everything works correctly
                @show n_nodes
                eigf[:,nfound] = tmp_eigf

                # update the number of solutions found
                nfound += 1
            end
            if high_derivative_workaround
                Estep = Estep*10
                high_derivative_workaround = false
            end
        end

        @label skip_interval
        prevdelta = delta
        E += Estep

        if (E > maximum(V))
           error("[Numerov] Could not find enough solutions with the parameters provided\n")
        end
    end

    return eigv, eigf
end


function findDelta!(E, V::Vector, centrifugal, h, xmin, xmax, bc_0_exp, bc_end_exp, strict_xc, verbose, yf, yb)

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
    xc = findEnergyIntersection(E,V,strict_xc)
    verbose && println("xc = ", xc)

    # calculation of the total potential term
    k2 = (E .- V) ./ h2m .- centrifugal

    numerov_forward!(h, xc, xmin, k2, yf)
    numerov_backward!(h, xc, xmax, k2, yb)

    delta = der5(yf,xc,h)/yf[xc] - der5(yb,xc,h)/yb[xc]

    verbose && @printf("yf[xc-1] = %0.9f,\tyf[xc] = %0.9f,\tyf[xc+1]=%0.9f\n", yf[xc-1], yf[xc], yf[xc+1])
    verbose && @printf("yb[xc-1] = %0.9f,\tyb[xc] = %0.9f,\tyb[xc+1]=%0.9f\n", yb[xc-1], yb[xc], yb[xc+1])
    verbose && @printf("derforward = %0.12f,  derback = %0.12f\n", der5(yf,xc,h)/yf[xc],der5(yb,xc,h)/yb[xc])
    verbose && @printf("E = %f\tDelta = %f\n", E, delta)

    return delta
end


# propagate the numerov solution up to xc+3
# xmin is used in order to not overwrite yf at the boundary provided to Numerov
@inbounds function numerov_forward!(h::Float64, xc::Int, xmin::Int, k2::Vector, yf::Vector)
    hh = h*h
    c0 = hh*k2[xmin]/12
    c_1 = hh*k2[xmin-1]/12
    c_2 = hh*k2[xmin-2]/12

    for x = xmin:xc+3
        yf[x] = (yf[x-1]*(2-10*c_1) - yf[x-2]*(1+c_2)) / (1+c0)
        c_2 = c_1
        c_1 = c0
        c0 = hh*k2[x+1]/12
    end
end

@inbounds function numerov_backward!(h::Float64, xc::Int, xmax::Int, k2::Vector, yb::Vector)
    hh = h*h
    c0 = hh*k2[xmax-2]/12
    c1 = hh*k2[xmax-1]/12
    c2 = hh*k2[xmax]/12

    for x = xmax-2:-1:xc-3
        yb[x] = (yb[x+1]*(2-10*c1) - yb[x+2]*(1+c2)) / (1+c0)
        #(x < xc+200) && println("yb = ", yb[x])
        c2 = c1
        c1 = c0
        c0 = hh*k2[x-1]/12
    end
end

function findEnergyIntersection(E::Float64, V::Vector, strict_xc::Bool)
    # _, xc = findmin(abs.(reverse(E .- V)))
    # xc = length(V) - xc
    # #xc, _ = secant(V .- E, (xmax*9)÷10, xmax-1, 10, 10^3)
    # #xc, _ = secant(V .- E, (xmax*9)÷10, xmax-1, 10, 10^3)
    # #xc, _ = secant(V .- E, 1, xmax-1, 10, 10^3)
    # if xc > length(V)*9÷10
    #     _, xc = findmin(abs.(E .- V))
    #     @warn "E-V intersection near the end of the domain, new one is at $xc"
    # end
    _, Vmin = findmin(V[1:length(V)*2÷3])
    _, xc = findmin(abs.(E .- V[Vmin+1:length(V)*3÷4]))
    xc += Vmin+1
    # if (xc<length(V)÷20)
    #     xc_old = xc
    #     _, xc = findmin(abs.(E .- V[xc_old+10:end]))
    #     xc += xc_old+10
    #     @warn "xc near 0 at $xc_old, choosing next intersection with potential in xc = $xc"
    # end
    if (xc<50)
        old_xc = xc
        _, Vmin = findmin(V[xc+50:length(V)*2÷3])
        Vmin += xc+50
        _, xc = findmin(abs.(E .- V[Vmin:length(V)*3÷4]))
        xc += Vmin
        @warn "xc was too small ($old_xc), found a new one at $xc"
    end

    if (strict_xc==true) && ((xc>length(V)*9÷10) || (xc<52))
        @warn "E-V intersection outside of the domain ($xc), E = $E, choosing xc in the middle of the domain.\nConsider checking the potential V for errors."
        xc = length(V)÷3
    end

    return xc
end

function count_nodes(X::Vector)
    nodes = 0
    numzeros = sum(X.==0)
    (sum(X.==0) > 0) && @info "the function has $numzeros points valued 0, check the number of nodes manually"
    @inbounds for i=1:length(X)-1
        if X[i]*X[i+1]<0
            nodes += 1
        end
    end
    return nodes
end


# TODO:
# return stuff even if E reached Vmax
# understand if it is ever useful to have longer than 2 boundary conditions, simplify code
# try again to use the secant method instead of the slow find() for the E-V intersection
# add test with a more complex potential
