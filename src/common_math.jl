
# finds the zero of f nearest to the interval [p1,p2]
function secant(f::Function, p1, p2, eps::Float64, nmax::Int)
    p = 0.
    for n = 1:nmax
        p = p2 - f(p2)*(p2-p1)/(f(p2)-f(p1))
        if abs(p-p2) < eps || f(p) == 0
            return p, f(p)
        end
        p1 = p2
        p2 = p
    end

    y = f(p)
    error("Method did not converge. The last iteration gives $p with function value $y")
end

# best used on convex functions, choose p1 and p2 carefully
function secant(f_::Vector, p1::Int, p2::Int, eps::Int, nmax::Int)
    p = 0.

    for n = 1:nmax # NOTE add inbounds when sure it works
        p = round(Int, p2 - f[p2]*(p2-p1)/(f[p2]-f[p1]))
        println(p)
        if (p<0 || p>length(f)) error("secant out of bounds at $p") end
        if abs(p-p2) < eps
            return p, f[p]
        end
        p1 = p2
        p2 = p
    end
    y = f[p]
    error("Method did not converge. The last iteration gives $p with function value $y")
end

function bettersecant(f_::Vector, p1::Int, p2::Int, eps::Int, nmax::Int)
    p = 0.
    # we create a symmetric version of f, to avoid it going below 0
    f = [reverse(f_); f_]
    p1 += length(f_)
    p2 += length(f_)
    for n = 1:nmax # NOTE add inbounds when sure it works
        println((f[p2]-f[p1]))
        println(f[p2])
        println(f[p2]*(p2-p1)/(f[p2]-f[p1]))
        println()
        p = round(Int, p2 - f[p2]*(p2-p1)/(f[p2]-f[p1]))
        println(p)
        if (p<0 || p>length(f)) error("secant out of bounds at $p") end
        if abs(p-p2) < eps
            return abs(p)-length(f_), f[abs(p)-length(f_)]
        end
        p1 = p2
        p2 = p
    end
    y = f[p]
    error("Method did not converge. The last iteration gives $p with function value $y")
end

# 3-point numerical derivative at point x
@inbounds function der3(F::Vector, x::Int, h::Float64)
    if x>1 && x < length(F)
        return (F[x+1]-F[x-1])/(2h)
    elseif x==1
        return (F[x+1]-F[x])/h
    else
        return (F[x]-F[x-1])/h
    end
end

# 5-point numerical derivative at point x
@inbounds function der5(F::Vector, x::Int, h::Float64)
    if x>2 && x < length(F)-1
        return (-F[x+2]+8*F[x+1]-8*F[x-1]+F[x-2])/(12h)
    else
        return der3(F,x,h)
    end
end

# 5-point numerical derivative of an array
@inbounds function der5(F::Vector, h::Float64)
    der = zeros(Float64,length(F))
    der[1] = (F[2]-F[1])/h
    der[2] = (F[3]-F[1])/(2h)
    for x = 3:length(F)-2
        der[x] = (-F[x+2]+8*F[x+1]-8*F[x-1]+F[x-2]) / (12h)
    end
    der[end-1] = (F[end]-F[end-2])/(2h)
    der[end] = (F[end]-F[end-1])/h

    return der
end



@inbounds function simpson_integral(f::Vector, xmax::Int, h)
    integral = h*(f[1]+f[xmax])/2

    for i = 2:2:xmax-1
        integral += h*(f[i-1] + 4*f[i] + f[i+1])/3
    end

    return integral
end

@inbounds function simpson_integral(f::Vector, h)
    integral = h*(f[1]+f[end])/2

    for i = 2:2:length(f)-1
        integral += h*(f[i-1] + 4*f[i] + f[i+1])/3
    end

    return integral
end
