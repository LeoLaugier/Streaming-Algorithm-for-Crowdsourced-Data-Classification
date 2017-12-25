function y = f(u,v,n)
    y = ( 1/(n-2) * sum(sqrt(v + 4*(n-1)/n^2 * (1-2*u))) )^2 ;