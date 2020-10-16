
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


# 5-point numerical derivative at point x
@inline der5(F::Vector, x::Int, h::Float64) = (-F[x+2]+8*F[x+1]-8*F[x-1]+F[x-2])/(12*h)


function simpson_integral(f::Vector, xmax::Int, h)
    integral = h*(f[1]+f[xmax])/2

    for i = 2:2:xmax-1
        integral += h*(f[i-1] + 4*f[i] + f[i+1])/3
    end

    return integral
end
