include("gauss_legendre.jl")
function t2ft_Q4(xnod_, lado, carga, espesor)

    #=Función que convierte las fuerzas superficiales aplicadas a un elemento
    finito rectangular de 4 nodos a sus correspondientes cargas nodales 
    equivalentes ft
    
    Recibe:
        xnod:  coordenadas nodales del elemento finito de 4 nodos
          xnod = [ x1e y1e
                   x2e y2e
                   x3e y3e
                   x4e y4e ]
        lado:  arista en la que se aplica la carga, puede tomar los siguientes
               valores: 12, 23, 34, 41
        carga: fuerza distribuida en los nodos
        
               [ t1x t1y t2x t2y ]; % si carga se aplica sobre lado 12
               [ t2x t2y t3x t3y ]; % si carga se aplica sobre lado 23
               [ t3x t3y t4x t4y ]; % si carga se aplica sobre lado 34
               [ t4x t4y t1x t1y ]; % si carga se aplica sobre lado 41
    
        espesor: espesor del elemento
    =#
    if     lado == 12   idx_ = [ 1 2 ]
    elseif lado == 23   idx_ = [ 2 3 ]
    elseif lado == 34   idx_ = [ 3 4 ]
    elseif lado == 41   idx_ = [ 4 5 ]
    else
        error("Únicamente se permiten los lados 12, 23, 34 o 41")
    end
    
    #xnod_ = xnod_[LaG[e,: ],:]
    nno = size(xnod_,1)
    if nno != 4
         error("Solo para elementos rectangulares de 4 nodos")
    end

     # parámetros para mejorar la lectura del código
     X, Y = 1, 2
   
     # se define el número de puntos de la cuadratura y se obtienen los puntos
     # de evaluación y los pesos de la cuadratura
     n_gl = 5                       # orden de la cuadratura de Gauss-Legendre
     x_gl, w_gl = gausslegendre_quad(n_gl)
     
     # se definen las funciones de forma unidimensionales y sus derivadas
     NN(xi)  = [ (1-xi)/2, (1+xi)/2 ]
     dNN_dxi(xi) = [     -1/2,      1/2 ]

    te_ = zeros(2*nno)
    te_[[1*idx_; 2*idx_]] = carga[:]

    ## Se calcula la integral
    suma   = zeros(2*nno,2*nno)
    N      = zeros(nno)
    dN_dxi_ = zeros(nno)

    for p = 1:n_gl
        N[idx_] = NN(x_gl[p])
        matN = zeros(2,2*nno)
        for i = 1:nno
            matN[:,[2*i-1, 2*i]] = [N[i] 0   
                                    0    N[i]]
        end
 
        dN_dxi_[idx_] = dNN_dxi(x_gl[p])
 
        dx_dxi_ = dot(dN_dxi_,xnod_[:,X])
        dy_dxi_ = dot(dN_dxi_,xnod_[:,Y])

        ds_dxi = hypot(dx_dxi_, dy_dxi_) 
 
        suma += matN'*matN*ds_dxi*w_gl[p]
        
     end
     ft = espesor.*suma*te_
    

    return  ft

end


function plot_def_esf_ang(xnod,esfdef, angulos, lab)

    X,Y = 1,2
 
    NL1, NL2, NL3, NL4= 1,2,3,4
    triangles = Vector{Vector{Int64}}(undef, 2*nef)
 
        # Para propósitos de graficación el EF se divide en 6 triángulos así: 
     #     
     #                             7 -------6--------5
     #                             |       /|\       |
     #                             | EFT6 / | \ EFT3 |
     #                             |     /  |  \     |
     #                             |    /   |   \    |
     #                             |   /    |    \   |
     #                             |  /     |     \  |
     #                             | /      |      \ |
     #                             8/  EFT5 | EFT4  \4
     #                             |\       |       /|
     #                             | \      |      / |
     #                             |  \     |     /  |
     #                             |   \    |    /   |
     #                             |    \   |   /    |
     #                             |     \  |  /     |
     #                             | EFT1 \ | / EFT2 |
     #                             |       \|/       |
     #                             1--------2--------3


    # Para propósitos de graficación con tripcolor, el EF MZC se divide en 2 triángulos así: 

    #                             4--------3
    #                             |\       |
    #                             | \ EFT2 |
    #                             |  \     |
    #                             |   \    |
    #                             |    \   |
    #                             |     \  |
    #                             | EFT1 \ |
    #                             |       \|
    #                             1--------2
       
 
     for e = 1:nef

        # se arma la matriz de correspondencia (LaG) de la malla
        triangles[2*e - 1] = LaG[e, [NL1, NL2, NL4]] .- 1
        triangles[2*e - 0] = LaG[e, [NL2, NL3, NL4]] .- 1
     end
 
       val_max = maximum(abs.(esfdef))
       fig, ax = subplots()
       # se grafica la malla de EFS, los colores en cada triángulo y las curvas 
       # de nivel
       for e = 1:nef
          # se dibujan las aristas
          nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                 plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.5, color = "gray")
       end
 
       im = ax.tripcolor(xnod[:, X], xnod[:, Y], triangles, esfdef,  cmap = "bwr",
                         shading = "gouraud", vmin = -val_max, vmax = val_max)
 
       ax.tricontour(xnod[:, X], xnod[:, Y], triangles, esfdef, 20)
       ax.tricontour(xnod[:, X], xnod[:, Y], triangles, esfdef, levels=[0], linewidths=3)
 
       fig.colorbar(im, ax = ax, format = "%6.3g")
 
       if ~isempty(angulos)
          # Grafique lineas que indican las direcciones principales de sigma_1
          norma = 1 # = esf si quiere proporcional
    
          for ang in angulos
             quiver(xnod[:,X],xnod[:,Y],              # En el nodo grafique una línea
                    norma.*cos.(ang),norma.*sin.(ang),# indicando la dirección
                    headlength=0,
                    headwidth = 0,
                    headaxislength = 0,
                    pivot="middle")
          end 
       end
       ylabel(lab)
       ax.set_aspect("equal")
       ax.autoscale(tight=true)
       tight_layout()
    return
 
end


