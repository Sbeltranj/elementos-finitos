
function t2ft_T3(xnod_EF, lado, carga, espesor)

    #=Convierte las fuerzas superficiales aplicadas a un EF triangular de 3
    nodos a sus correspondientes cargas nodales equivalentes ft
    xnod_EF =          ([[x1e, y1e]
                        [x2e, y2e]
                        [x3e, y3e]])
    lado = 12, 23, 31
    carga = [ t1x, t1y, t2x, t2y ]   # si la carga se aplica sobre el lado 12
            [ t2x, t2y, t3x, t3y ]   # si la carga se aplica sobre el lado 23
            [ t3x, t3y, t1x, t1y ]   # si la carga se aplica sobre el lado 31
    =#

    if lado == 12
        # Fuerzas sobre el lado 12
        # Se calcula la longitud del lado 12
        L12 = hypot(xnod_EF[NL1,X] - xnod_EF[NL2,X], xnod_EF[NL1,Y] - xnod_EF[NL2,Y])

        # Fuerzas distribuidas aplicadas en los nodos 1 y 2 locales
        t1x, t1y, t2x, t2y = carga
        ft =            espesor*([ (L12*(2*t1x + t2x))/6
                                   (L12*(2*t1y + t2y))/6
                                   (L12*(t1x + 2*t2x))/6
                                   (L12*(t1y + 2*t2y))/6
                                    0
                                    0                      ])
    elseif lado == 23
        # Fuerzas sobre el lado 23
        # Se calcula la longitud del lado 23
        L23 =  hypot(xnod_EF[NL2,X] - xnod_EF[NL3,X], xnod_EF[NL2,Y] - xnod_EF[NL3,Y])

        # Fuerzas distribuidas aplicadas en los nodos 2 y 3 locales
        t2x, t2y, t3x, t3y = carga
        ft =           espesor*([ 0
                                  0
                                (L23*(2*t2x + t3x))/6
                                (L23*(2*t2y + t3y))/6
                                (L23*(t2x + 2*t3x))/6
                                (L23*(t2y + 2*t3y))/6 ])
    elseif lado == 31
        # Fuerzas sobre el lado 31
        # Se calcula la longitud del lado 31
        L31 = hypot(xnod_EF[NL3,X] - xnod_EF[NL1,X], xnod_EF[NL3,Y] - xnod_EF[NL1,Y])

        # Fuerzas distribuidas aplicadas en los nodos 3 y 1 locales
        t3x, t3y, t1x, t1y = carga
        ft =         espesor*([ (L31*(2*t1x + t3x))/6
                                (L31*(2*t1y + t3y))/6
                                0
                                0
                                (L31*(t1x + 2*t3x))/6
                                (L31*(t1y + 2*t3y))/6 ])
    else
        println("error")
    end


    return ft

end
#=
# NOTA: Los vectores se calcularon con el siguiente programa de MATLAB:
# Elemento triangular de 3 nodos
clear
clc
syms N1 N2 N3
syms t1x t1y t2x t2y t3x t3y
syms xi L12 L23 L31
N = [ N1 0   N2 0   N3 0
  0  N1  0  N2  0  N3 ];
t = [ t1x t1y t2x t2y t3x t3y ].';
res = N.'*N*t;
res12 = res;
res12 = subs(res12, N1, (1-xi)/2);
res12 = subs(res12, N2, (1+xi)/2);
res12 = subs(res12, N3, 0);
ds_dxi = L12/2;
f12 = int(res12*ds_dxi, xi, -1 ,1)
res23 = res;
res23 = subs(res23, N2, (1-xi)/2);
res23 = subs(res23, N3, (1+xi)/2);
res23 = subs(res23, N1, 0);
ds_dxi = L23/2;
f23 = int(res23*ds_dxi, xi, -1 ,1)
res31 = res;
res31 = subs(res31, N3, (1-xi)/2);
res31 = subs(res31, N1, (1+xi)/2);
res31 = subs(res31, N2, 0);
ds_dxi = L31/2;
f31 = int(res31*ds_dxi, xi, -1 ,1)
'''
=#

function plot_def_esf_ang(xnod,esfdef, angulos, lab)

    #Actualizando a tripcolor:
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/tripcolor_demo.html
       fig, ax = subplots()
 
       X, Y = 1,2
       NL1, NL2, NL3 = 1, 2, 3
 
       triangles = Vector{Vector{Int64}}(undef, nef)
 
       for e = 1:nef
 
         # se arma la matriz de correspondencia (LaG) de la nueva malla triangular
         triangles[e] = LaG[e, [NL1, NL2, NL3]] .- 1
   
       end

       #Visitar :https://github.com/diegoandresalvarez/elementosfinitos/blob/master/diapositivas/05a_EF_2D_T3_y_Q4.pdf
 
 
     #                                     2
     #                                    / \       
     #                                   /   \      
     #                                  /     \     
     #                                 /       \   
     #                                /         \   
     #                               /           \  
     #                              /             \ 
     #                            3/  ---- | ----  \1
 
 
         # se grafica la malla de EFS, los colores en cada triángulo y las curvas
         # de nivel
 
         ax.triplot(xnod[:,X], xnod[:,Y], triangles,  lw=0.5, color="gray")
 
         # se promedian los esfuerzos y las deformaciones en los nodos de modo
         # que en la variable "var" se encuentren los valores alisados del
         # esfuerzo o de la deformación a graficar
         
         num_elem_ady = zeros(nno)
 
 
         for e = 1:nef
             num_elem_ady[LaG[e,:]] .+= 1
         end
 
         var = zeros(nno)
         for e = 1:nef
             var[LaG[e,:]] .+= esfdef[e]
         end
         
     
         var  = var./num_elem_ady
 
         # se encuentra el máximo en valor absoluto para ajustar el colorbar()
         val_max = maximum(var)

         
 
         im = ax.tripcolor(xnod[:, X], xnod[:, Y], triangles, var,  cmap = "bwr",
                          shading = "gouraud", vmin = -val_max, vmax = val_max)
 
         fig.colorbar(im, ax = ax, format = "%6.3g")
         ax.tricontour(xnod[:, X], xnod[:, Y], triangles, var, 20)
         axis("equal")
 
         ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
       if ~isempty(angulos)
          # Grafique lineas que indican las direcciones principales de sigma_1
          norma = 1 
    
          for ang in angulos
             quiver(cgx[:,X],cgy[:,X],                 # En el nodo grafique una línea
                    norma.*cos.(ang), norma.*sin.(ang),# indicando la dirección
                    headlength=0,
                    headwidth = 0,
                    headaxislength = 0,
                    pivot="middle")
          end
          # scatter(xnod[:,X],xnod[:,Y], color="k", s=1)  # para poner el punto
          # http://matplotlib.org/examples/pylab_examples/quiver_demo.html
       end
       ylabel(lab)
       ax.set_aspect("equal")
       ax.autoscale(tight=true)
    return
 
end
#%% Fin 