function y = h(u,v,n)
    b = zeros(n,1);
    for i = 1:n
        b(i) = 1/2 + n/4 * (sqrt(v + 4*(n-1)/n^2 * (1-2*u(i))) - sqrt(v)) ;
    end
    y = b ;