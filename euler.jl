## Euler integrator
function euler(f,y0,t)
    d = prod(size(y0))
    n = length(t)

    sol = zeros(n,d)
    sol[1,:] = y0
    
    for k=2:n
        dt = t[k]-t[k-1]
        sol[k,:] = sol[k-1,:] + dt*f(sol[k-1,:],t[k-1])
    end
    return sol
end
