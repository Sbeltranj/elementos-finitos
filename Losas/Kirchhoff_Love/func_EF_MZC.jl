##Programa elaborado en JULIA 1.6.3

#Santiago Beltrán Jaramillo
#sbeltran@unal.edu.co

#Función para dibujar los diferentes gráficos del EF MZC,
#para la losa de Kirchhoff_Love
#Tener instalada la librería de matplotlib (PYTHON-pip install matplotlib)
#además Pkg.add("PyPlot") en consola de JULIA

using PyPlot

function plot_def_esf_ang(xnod,esfdef, lab)
    
    #Constantes para facilitar el código
    X,Y = 1,2
    
    # EF MZC de 4 nodos.
    NL1, NL2, NL3, NL4 = 1,2,3,4

    #separó memoria para el triángulo de graficación.
    triangles = Vector{Vector{Int64}}(undef, 2*nef)

    for e = 1:nef

        # se arma la matriz de correspondencia (LaG) de la malla
        triangles[2*e - 1] = LaG[e, [NL1, NL2, NL4]] .- 1
        triangles[2*e - 0] = LaG[e, [NL2, NL3, NL4]] .- 1
     end
 
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
       
      
     fig, ax = subplots(figsize=(10, 10))

     for e = 1:length(esfdef)
            
            if e == 1

                #creo un subplot de 2f, 3col.
                subplot(231)

                #se hace la interpolación en los triángulos.
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, esfdef[e],  cmap = "bwr",
                          shading = "gouraud")
                
                #algunos ajustes a la gráfica.
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar() #barra de colores
                title(lab[e])

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                 end

            elseif e == 2

                #creo un subplot de 2f, 3col.
                subplot(232)

                #se hace la interpolación en los triángulos.
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, esfdef[e],  cmap = "bwr",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar()
                title(lab[e])

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                 end

            elseif e == 3

                subplot(233)
                #se hace la interpolación en los triángulos.
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, esfdef[e],  cmap = "bwr",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar()
                title(lab[e])

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                 end
            else 
                error("                                      ")
                
            end

       end

    return
 
end