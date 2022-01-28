function Bs_RM(xi, eta, xe, ye, Nforma, dN_dxi, dN_deta)
    ## Calcula la matriz de deformaciones por cortante Bs
    #
    # Bs_RM(xi, eta, xe, ye, Nforma, dN_dxi, dN_deta)
    # 
    # (xi, eta)        punto de integración de GL
    # (xe, ye)         coordenadas de los nodos del EF
    # Nforma           function  con funciones de forma del EF
    # dN_dxi, dN_deta  function  derivadas de las funciones de forma del EF
    
    ## Se evalúan las funciones de forma en los puntos de integración
    ## de Gauss-Legendre
    NNforma = Nforma(xi, eta);
    
    ## Se evalúan las derivadas de las funciones de forma en los puntos
    ## de integracion de Gauss-Legendre
    ddN_dxi  = dN_dxi(xi, eta);
    ddN_deta = dN_deta(xi, eta);
    
    ## Se utilizan las funciones de forma de w para el calculo de la 
    ## transformación isoparamétrica
    dx_dxi  = sum(ddN_dxi  .* xe);   dy_dxi  = sum(ddN_dxi  .* ye);
    dx_deta = sum(ddN_deta .* xe);   dy_deta = sum(ddN_deta .* ye);
    
    ## Se ensambla la matriz Jacobiana del elemento
    Je = [ dx_dxi   dy_dxi
           dx_deta  dy_deta ];
    
    ## Se calcula el determinante del Jacobiano
    det_Je = det(Je);
    if det_Je <= 0
       error("El det_Je es negativo");
    end
    
    ## Se ensambla la matriz de deformacion del elemento Bs
    nno = length(xe);
    Bs  = zeros(2,3*nno);
    for i = 1:nno
       dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je;
       dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je;
    
       # se ensambla y asigna la matriz Bs_i
       Bs[:,(3*i-2):(3*i)] = [ dNi_dx  -NNforma[i]  0     
                               dNi_dy  0            -NNforma[i] ];                        
    end
    
    return Bs, det_Je, Je

end