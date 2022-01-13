include("gausslegendre_quad.jl")

function gausslegendre_quad_hexa(n_gl)

    x, w = gausslegendre_quad(n_gl);
    
    x_gl = zeros(n_gl^3,3);
    w_gl = zeros(n_gl^3,1);
    
    r = 0;
    for i = 1:n_gl
       for j = 1:n_gl
          for k = 1:n_gl
             r = r+1;
             x_gl[r,:] = [ x[i] x[j] x[k] ];
             w_gl[r] = w[i]*w[j]*w[k];
          end
       end
    end
    
    return x_gl, w_gl
end

function N_H20(xi,eta,zeta)
    # Funciones de forma del EF hexaédrico serendípito isoparamétrico de 20 nodos
    
    N = [ 
         ((eta - 1)*(xi - 1)*(zeta - 1)*(eta + xi + zeta + 2))/8     # N1
        -((xi^2 - 1)*(eta - 1)*(zeta - 1))/4                         # N2
        -((eta - 1)*(xi + 1)*(zeta - 1)*(eta - xi + zeta + 2))/8     # N3
         ((eta^2 - 1)*(xi + 1)*(zeta - 1))/4                         # N4
        -((eta + 1)*(xi + 1)*(zeta - 1)*(eta + xi - zeta - 2))/8     # N5
         ((xi^2 - 1)*(eta + 1)*(zeta - 1))/4                         # N6
        -((eta + 1)*(xi - 1)*(zeta - 1)*(xi - eta + zeta + 2))/8     # N7
        -((eta^2 - 1)*(xi - 1)*(zeta - 1))/4                         # N8
        -((zeta^2 - 1)*(eta - 1)*(xi - 1))/4                         # N9
         ((zeta^2 - 1)*(eta - 1)*(xi + 1))/4                         # N10
        -((zeta^2 - 1)*(eta + 1)*(xi + 1))/4                         # N11
         ((zeta^2 - 1)*(eta + 1)*(xi - 1))/4                         # N12
        -((eta - 1)*(xi - 1)*(zeta + 1)*(eta + xi - zeta + 2))/8     # N13
         ((xi^2 - 1)*(eta - 1)*(zeta + 1))/4                         # N14
         ((eta - 1)*(xi + 1)*(zeta + 1)*(eta - xi - zeta + 2))/8     # N15
        -((eta^2 - 1)*(xi + 1)*(zeta + 1))/4                         # N16
         ((eta + 1)*(xi + 1)*(zeta + 1)*(eta + xi + zeta - 2))/8     # N17
        -((xi^2 - 1)*(eta + 1)*(zeta + 1))/4                         # N18
        -((eta + 1)*(xi - 1)*(zeta + 1)*(eta - xi + zeta - 2))/8     # N19
         ((eta^2 - 1)*(xi - 1)*(zeta + 1))/4                         # N20
    ];                                                              
    
    return N
end

function dN_dxi_H20(xi,eta,zeta)
    # Derivadas de las funciones de forma con respecto a la variable "xi" del
    # EF hexaédrico serendípito isoparamétrico de 20 nodos
    
    dN_dxi = [ 
      ((eta - 1)*(zeta - 1)*(eta + 2*xi + zeta + 1))/8   # dN1_dxi 
     -(xi*(eta - 1)*(zeta - 1))/2                        # dN2_dxi 
     -((eta - 1)*(zeta - 1)*(eta - 2*xi + zeta + 1))/8   # dN3_dxi 
      ((eta^2 - 1)*(zeta - 1))/4                         # dN4_dxi 
     -((eta + 1)*(zeta - 1)*(eta + 2*xi - zeta - 1))/8   # dN5_dxi 
      (xi*(eta + 1)*(zeta - 1))/2                        # dN6_dxi 
     -((eta + 1)*(zeta - 1)*(2*xi - eta + zeta + 1))/8   # dN7_dxi 
     -((eta^2 - 1)*(zeta - 1))/4                         # dN8_dxi 
     -((zeta^2 - 1)*(eta - 1))/4                         # dN9_dxi 
      ((zeta^2 - 1)*(eta - 1))/4                         # dN10_dxi
     -((zeta^2 - 1)*(eta + 1))/4                         # dN11_dxi
      ((zeta^2 - 1)*(eta + 1))/4                         # dN12_dxi
     -((eta - 1)*(zeta + 1)*(eta + 2*xi - zeta + 1))/8   # dN13_dxi
      (xi*(eta - 1)*(zeta + 1))/2                        # dN14_dxi
      ((eta - 1)*(zeta + 1)*(eta - 2*xi - zeta + 1))/8   # dN15_dxi
     -((eta^2 - 1)*(zeta + 1))/4                         # dN16_dxi
      ((eta + 1)*(zeta + 1)*(eta + 2*xi + zeta - 1))/8   # dN17_dxi
     -(xi*(eta + 1)*(zeta + 1))/2                        # dN18_dxi
     -((eta + 1)*(zeta + 1)*(eta - 2*xi + zeta - 1))/8   # dN19_dxi
      ((eta^2 - 1)*(zeta + 1))/4                         # dN20_dxi
    ];
    
    return dN_dxi 

end


function dN_deta_H20(xi,eta,zeta)
    # Derivadas de las funciones de forma con respecto a la variable "eta" del
    # EF hexaédrico serendípito isoparamétrico de 20 nodos
    
    dN_deta = [ 
         ((xi - 1)*(zeta - 1)*(2*eta + xi + zeta + 1))/8    # dN1_deta  
        -((xi^2 - 1)*(zeta - 1))/4                          # dN2_deta  
        -((xi + 1)*(zeta - 1)*(2*eta - xi + zeta + 1))/8    # dN3_deta  
         (eta*(xi + 1)*(zeta - 1))/2                        # dN4_deta  
        -((xi + 1)*(zeta - 1)*(2*eta + xi - zeta - 1))/8    # dN5_deta  
         ((xi^2 - 1)*(zeta - 1))/4                          # dN6_deta  
        -((xi - 1)*(zeta - 1)*(xi - 2*eta + zeta + 1))/8    # dN7_deta  
        -(eta*(xi - 1)*(zeta - 1))/2                        # dN8_deta  
        -((zeta^2 - 1)*(xi - 1))/4                          # dN9_deta  
         ((zeta^2 - 1)*(xi + 1))/4                          # dN10_deta 
        -((zeta^2 - 1)*(xi + 1))/4                          # dN11_deta 
         ((zeta^2 - 1)*(xi - 1))/4                          # dN12_deta 
        -((xi - 1)*(zeta + 1)*(2*eta + xi - zeta + 1))/8    # dN13_deta 
         ((xi^2 - 1)*(zeta + 1))/4                          # dN14_deta 
         ((xi + 1)*(zeta + 1)*(2*eta - xi - zeta + 1))/8    # dN15_deta 
        -(eta*(xi + 1)*(zeta + 1))/2                        # dN16_deta 
         ((xi + 1)*(zeta + 1)*(2*eta + xi + zeta - 1))/8    # dN17_deta 
        -((xi^2 - 1)*(zeta + 1))/4                          # dN18_deta 
        -((xi - 1)*(zeta + 1)*(2*eta - xi + zeta - 1))/8    # dN19_deta 
         (eta*(xi - 1)*(zeta + 1))/2                        # dN20_deta 
    ]
    
    return dN_deta 

end

function dN_dzeta_H20(xi,eta,zeta)
    # Derivadas de las funciones de forma con respecto a la variable "zeta" del
    # EF hexaedrico serendipito isoparametrico de 20 nodos
    
    dN_dzeta = [ 
      ((eta - 1)*(xi - 1)*(eta + xi + 2*zeta + 1))/8    # dN1_dzeta 
     -((xi^2 - 1)*(eta - 1))/4                          # dN2_dzeta 
     -((eta - 1)*(xi + 1)*(eta - xi + 2*zeta + 1))/8    # dN3_dzeta 
      ((eta^2 - 1)*(xi + 1))/4                          # dN4_dzeta 
     -((eta + 1)*(xi + 1)*(eta + xi - 2*zeta - 1))/8    # dN5_dzeta 
      ((xi^2 - 1)*(eta + 1))/4                          # dN6_dzeta 
     -((eta + 1)*(xi - 1)*(xi - eta + 2*zeta + 1))/8    # dN7_dzeta 
     -((eta^2 - 1)*(xi - 1))/4                          # dN8_dzeta 
     -(zeta*(eta - 1)*(xi - 1))/2                       # dN9_dzeta 
      (zeta*(eta - 1)*(xi + 1))/2                       # dN10_dzeta
     -(zeta*(eta + 1)*(xi + 1))/2                       # dN11_dzeta
      (zeta*(eta + 1)*(xi - 1))/2                       # dN12_dzeta
     -((eta - 1)*(xi - 1)*(eta + xi - 2*zeta + 1))/8    # dN13_dzeta
      ((xi^2 - 1)*(eta - 1))/4                          # dN14_dzeta
      ((eta - 1)*(xi + 1)*(eta - xi - 2*zeta + 1))/8    # dN15_dzeta
     -((eta^2 - 1)*(xi + 1))/4                          # dN16_dzeta
      ((eta + 1)*(xi + 1)*(eta + xi + 2*zeta - 1))/8    # dN17_dzeta
     -((xi^2 - 1)*(eta + 1))/4                          # dN18_dzeta
     -((eta + 1)*(xi - 1)*(eta - xi + 2*zeta - 1))/8    # dN19_dzeta
      ((eta^2 - 1)*(xi - 1))/4                          # dN20_dzeta
    ];
    
    return dN_dzeta 
end

function matriz_extrapolacion_esfuerzos_H20()

    A = [ 
(3*3^(1/2))/4 + 5/4      -3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4
    3^(1/2)/4 + 1/2                   -1/4                   -1/4        1/2 - 3^(1/2)/4        3^(1/2)/4 + 1/2                   -1/4                   -1/4        1/2 - 3^(1/2)/4
  - 3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4    (3*3^(1/2))/4 + 5/4      -3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4
               -1/4        1/2 - 3^(1/2)/4                   -1/4        1/2 - 3^(1/2)/4        3^(1/2)/4 + 1/2                   -1/4        3^(1/2)/4 + 1/2                   -1/4
    3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4      -3^(1/2)/4 - 1/4
               -1/4        1/2 - 3^(1/2)/4        3^(1/2)/4 + 1/2                   -1/4                   -1/4        1/2 - 3^(1/2)/4        3^(1/2)/4 + 1/2                   -1/4
  - 3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4
    3^(1/2)/4 + 1/2                   -1/4        3^(1/2)/4 + 1/2                   -1/4                   -1/4        1/2 - 3^(1/2)/4                   -1/4        1/2 - 3^(1/2)/4
    3^(1/2)/4 + 1/2        3^(1/2)/4 + 1/2                   -1/4                   -1/4                   -1/4                   -1/4        1/2 - 3^(1/2)/4        1/2 - 3^(1/2)/4
               -1/4                   -1/4        1/2 - 3^(1/2)/4        1/2 - 3^(1/2)/4        3^(1/2)/4 + 1/2        3^(1/2)/4 + 1/2                   -1/4                   -1/4
    1/2 - 3^(1/2)/4        1/2 - 3^(1/2)/4                   -1/4                   -1/4                   -1/4                   -1/4        3^(1/2)/4 + 1/2        3^(1/2)/4 + 1/2
               -1/4                   -1/4        3^(1/2)/4 + 1/2        3^(1/2)/4 + 1/2        1/2 - 3^(1/2)/4        1/2 - 3^(1/2)/4                   -1/4                   -1/4
  - 3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4        3^(1/2)/4 - 1/4
               -1/4        3^(1/2)/4 + 1/2        1/2 - 3^(1/2)/4                   -1/4                   -1/4        3^(1/2)/4 + 1/2        1/2 - 3^(1/2)/4                   -1/4
    3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4    5/4 - (3*3^(1/2))/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4
    1/2 - 3^(1/2)/4                   -1/4        1/2 - 3^(1/2)/4                   -1/4                   -1/4        3^(1/2)/4 + 1/2                   -1/4        3^(1/2)/4 + 1/2
5/4 - (3*3^(1/2))/4        3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4
    1/2 - 3^(1/2)/4                   -1/4                   -1/4        3^(1/2)/4 + 1/2        1/2 - 3^(1/2)/4                   -1/4                   -1/4        3^(1/2)/4 + 1/2
    3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4    (3*3^(1/2))/4 + 5/4    5/4 - (3*3^(1/2))/4        3^(1/2)/4 - 1/4        3^(1/2)/4 - 1/4      -3^(1/2)/4 - 1/4
               -1/4        3^(1/2)/4 + 1/2                   -1/4        3^(1/2)/4 + 1/2        1/2 - 3^(1/2)/4                   -1/4        1/2 - 3^(1/2)/4                   -1/4 ];
  return A
end