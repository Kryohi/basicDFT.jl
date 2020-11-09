
# finds the zero of f nearest to the interval [p1,p2]
function secant(f::Function, p1, p2, eps::Float64, nmax::Int; MT=false)
    p = 0.
    y1, y2, y = 0., 0., 0.
    for n = 1:nmax
        if MT
            pp = [p1;p2]
            yy = [0.;0.]
            Threads.@threads for i=1:2
                yy[i] = f(pp[i])
            end
            y1, y2 = yy[1], yy[2]
        else
            y1 = f(p1)
            y2 = f(p2)
        end

        if isnan(y1) || isnan(y2)
            @warn "Secant returning NaNs, with p1=$p1, p2=$p2"
            return NaN, NaN
        end

        p = p2 - y2*(p2-p1)/(y2-y1)

        if abs(p-p2) < eps || y2 == 0
            y = f(p)
            isnan(y) && @warn "Secant returning NaNs, with p=$p, p2=$p2"
            return p, y
        end
        p1 = p2
        p2 = p
    end
    @warn "Method did not converge. The last iteration gives $p with function value $y"
    return NaN, NaN
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
    @warn "Method did not converge. The last iteration gives $p with function value $y"
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
    @warn "Method did not converge. The last iteration gives $p with function value $y"
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
    if x>2 && x<length(F)-1
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
    integral = 0.

    for i = 2:2:xmax-1
        integral += h*(f[i-1] + 4*f[i] + f[i+1])/3
    end

    return integral
end

@inbounds function simpson_integral(f::Vector, h)
    integral = 0.

    for i = 2:2:length(f)-1
        integral += h*(f[i-1] + 4*f[i] + f[i+1])/3
    end

    return integral
end


## filters & interpolation (using DSP.jl would be likely better)
# could be much faster but it's not important
function moving_average(X::Vector, w::Float64, h::Float64)
      N = length(X)
      Xsmooth = zeros(N)
      l = ceil(Int,w/h)
      Xsmooth[1:l] = X[1:l]
      Xsmooth[N-l:N] = X[N-l:N]
      for i=l+1:N-l
            Xsmooth[i] = mean(X[i-l:i+l])
      end
      return Xsmooth
end
function exp_moving_average(X::Vector, w::Float64, h::Float64)
      N = length(X)
      Xsmooth = zeros(N)
      l = ceil(Int,w/h)
      Xsmooth[1:l] = X[1:l]
      Xsmooth[N-l:N] = X[N-l:N]
      tail = exp.(-Vector(1:l) .* (1/4))
      c = [reverse(tail); 1.; tail]
      @show c /= sum(c)
      for i=l+1:N-l-1
            for j=i-l:i+l
                  Xsmooth[i] += c[j-i+l+1]*X[j]
            end
      end
      return Xsmooth
end
# fixes oscillations of period 2h
function custom_moving_average(X::Vector)
    N = length(X)
    Xsmooth = zeros(N)
    Xsmooth[1] = X[1]
    Xsmooth[N] = X[N]
    for i=2:N-1
      Xsmooth[i] = 0.5*X[i] + 0.25*X[i-1] + 0.25*X[i+1]
    end
    return Xsmooth
end
