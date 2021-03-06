function tetra_quad(n)
    # Cuadratura de Gauss-Legendre para un tetraedro
       if n == 1 # Precision lineal
          L2L3L4 = [0.25 0.25 0.25];
          W = 1;
       elseif n == 4 # Precision cuadrática
          a = 0.58541020; b = 0.13819660;
          L2L3L4 = [ a b b
                     b a b
                     b b a
                     b b b ];
          W = [0.25; 0.25; 0.25; 0.25];
       elseif n == 5 # Precisión cúbica
          L2L3L4 = [0.25 0.25 0.25
                    1/3  1/6  1/6
                    1/6  1/3  1/6
                    1/6  1/6  1/3
                    1/6  1/6  1/6 ];      
          W = [-0.8; 0.45; 0.45; 0.45; 0.45];
       else
          error("Numero de puntos de la cuadratura no soportada");         
       end
       
     W = W/6;
    
    return L2L3L4,W
end