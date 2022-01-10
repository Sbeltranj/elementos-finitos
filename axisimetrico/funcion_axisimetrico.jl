include("gausslegendre_quad.jl")

function t2ft_R89(xnod_, lado, carga)
    # Esta función convierte las fuerzas superficiales aplicadas a un elemento
    # finito rectangular de 8 (serendipito) o 9 (lagrangiano) nodos a sus
    # correspondientes cargas nodales equivalentes ft
    #
    # SERENDIPITO 8          LAGRANGIANO 9
    # xnod = [ x1e y1e       xnod = [ x1e y1e
    #          x2e y2e                x2e y2e
    #          ... ...                ... ...
    #          x8e y8e ];             x9e y9e ];
    #
    # lado = 123, 345, 567, 781
    #
    # carga = [ t1x t1y t2x t2y t3x t3y ]; # si carga se aplica sobre lado 123
    #         [ t3x t3y t4x t4y t5x t5y ]; # si carga se aplica sobre lado 345
    #         [ t5x t5y t6x t6y t7x t7y ]; # si carga se aplica sobre lado 567
    #         [ t7x t7y t8x t8y t1x t1y ]; # si carga se aplica sobre lado 781
 
    ## Se definen algunas constantes
    X = 1; Y = 2

    nno = size(xnod_,1)

    if nno != 8
         error("Solo para elementos rectangulares de 8 nodos")
    end
 
    ## párametros de la cuadratura de Gauss-Legendre
    n_gl = 2                       # orden de la cuadratura de Gauss-Legendre
    x_gl, w_gl = gausslegendre_quad(n_gl)
 
    ## Se definen las funciones de forma unidimensionales y sus derivadas
    NN(xi) = [
                xi .*(xi .-1)/2       # N1
                (1 .+ xi).*(1 .-xi)   # N2
                xi .*(1 .+xi)/2    ]     # N3
 
    dNN_dxi(xi) = [
                    xi .- 1/2           # dN1_dxi
                    -2 .* xi             # dN2_dxi
                    xi .+ 1/2 ]         # dN3_dxi
 
    ## Se definen los indices de los lados
    if     lado == 123   idx = [ 1 2 3 ]
    elseif lado == 345   idx = [ 3 4 5 ]
    elseif lado == 567   idx = [ 5 6 7 ]
    elseif lado == 781   idx = [ 7 8 1 ]
    else
    #error("Unicamente se permiten los lados 123, 345, 567 o 781")
    end
 
    ## Se calcula el vector de fuerzas distribuidas en los nodos
    te = zeros(2*nno)
    te[[2*idx.-1; 2*idx]] = carga[:]

    ## Se calcula la integral
    suma    = zeros(2*nno,2*nno)
    N       = zeros(nno)
    dN_dxi_ = zeros(nno)

    for p = 1:n_gl
        N[idx] = NN(x_gl[p])
        matN = zeros(2,2*nno)
        for i = 1:nno
            matN[:,[2*i-1, 2*i]] = [N[i] 0   
                                    0    N[i]]
        end
 
        dN_dxi_[idx] = dNN_dxi(x_gl[p])
              r       = dot(N,      xnod_[:,X])
              dx_dxi_ = dot(dN_dxi_,xnod_[:,X])
              dy_dxi_ = dot(dN_dxi_,xnod_[:,Y])
              ds_dxi  = hypot(dx_dxi_, dy_dxi_) 

        # y se calcula la sumatoria

        suma += r * matN'*matN*ds_dxi*w_gl[p]
    end
 
    ft = suma*te
 
    return ft
 end


function plot_def_esf_ang(xnod,esfdef, angulos, lab)

   X,Y = 1,2

   NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = 1,2,3,4,5,6,7,8
   triangle = Vector{Vector{Int64}}(undef, 6*nef)

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

   for e = 1:nef

      # se arma la matriz de correspondencia (LaG) de la nueva malla triangular
      triangle[6*e - 5] = LaG[e, [NL1, NL2, NL8]] .- 1
      triangle[6*e - 4] = LaG[e, [NL2, NL3, NL4]] .- 1
      triangle[6*e - 3] = LaG[e, [NL4, NL5, NL6]] .- 1
      triangle[6*e - 2] = LaG[e, [NL2, NL4, NL6]] .- 1
      triangle[6*e - 1] = LaG[e, [NL2, NL6, NL8]] .- 1
      triangle[6*e - 0] = LaG[e, [NL6, NL7, NL8]] .- 1

   end

      val_max = maximum(abs.(esfdef))
      fig, ax = subplots()
       # se grafica la malla de EFS, los colores en cada triángulo y las curvas 
       # de nivel
      for e = 1:nef
         # se dibujan las aristas
         nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
                plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.5, color = "gray")
      end

      im = ax.tripcolor(xnod[:, X], xnod[:, Y], triangle, esfdef,  cmap = "bwr",
                        shading = "gouraud", vmin = -val_max, vmax = val_max)

      ax.tricontour(xnod[:, X], xnod[:, Y], triangle, esfdef, 20)
      ax.tricontour(xnod[:, X], xnod[:, Y], triangle, esfdef, levels=[0], linewidths=3)

      fig.colorbar(im, ax = ax, format = "%6.3g")

      if ~isempty(angulos)
         # Grafique lineas que indican las direcciones principales de sigma_1
         norma = 1 # = esf si quiere proporcional
   
         for ang in angulos
            quiver(xnod[:,X],xnod[:,Y],              # En el nodo grafique una línea
                   esfdef.*cos.(ang),esfdef.*sin.(ang),# indicando la dirección
                   headlength=0,
                   headwidth = 0,
                   headaxislength = 0,
                   pivot="middle")
         end
      end
     # ylabel(lab)
      ax.set_xlabel("r [m]")
      ax.set_ylabel("z [m]")
      ax.set_aspect("equal")
      ax.autoscale(tight=true)
      ax.set_title(lab, fontsize=20)
      tight_layout()
   return

end


