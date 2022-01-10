using Polynomials

function gausslegendre_quad(m)
    # Integration using a Gauss-Legendre quadrature
    
    ## Calculation of the Legendre polynomials using Bonnet's recursion:
    #           P_n(x) = ((2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x))/n
    
    # Remember that  JULIA does not make 0-based indexing of arrays
    
    P = Vector{Polynomial{Float64}}(undef, m+1)
    P[0 + 1] =     Polynomial([1.0])       # P_{0}(x) = 1
    P[1 + 1] = x = Polynomial([0.0, 1.0])  # P_{1}(x) = x
    
    for n = 2:(m-1 + 1)
        P[n + 1] = ((2*n - 1)*x*P[n-1 + 1] - (n-1)*P[n-2 + 1])/n
    end
    
        ## Roots
        xi = sort(roots(Polynomial(P[m+1])));
    
        ## Weights
        s = derivative(Polynomial(P[m+1]));
    
        w = 2.0 ./ ((1 .- xi.^2).*(s.(xi)).^2);
    
    return xi, w
end

