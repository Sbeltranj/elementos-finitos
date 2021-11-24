##Programa elaborado en JULIA 1.6.3

#Santiago Beltrán Jaramillo
#sbeltran@unal.edu.co

#Función para dibujar los diferentes gráficos del EF MZC,
#para la losa de Kirchhoff_Love
#Tener instalada la librería de matplotlib (PYTHON-pip install matplotlib)
#además Pkg.add("PyPlot") en consola de JULIA

function plot_mom_Q_ang(xnod, mom_Q, ang_mf1, ang_mf2, ang_mt, lab)

    X,Y = 1,2
 
    NL1, NL2, NL3, NL4 = 1,2,3,4
    
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

       for e = 1:length(mom_Q)
            
            if e == 1
                subplot(231)
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "jet",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar()
                title(lab[e])

                if ~isempty(ang_mf1)

                    # Grafique lineas que indican las direcciones principales
                    norma = 1

                    for ang1 in ang_mf1
                        quiver(xnod[:,X],xnod[:,Y],                   # En el nodo grafique una línea
                                norma.*cos.(ang1), norma.*sin.(ang1), # indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end
               
                # se dibujan las aristas
                for e = 1:nef
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                 end

            elseif e == 2
                subplot(232)
                
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "jet",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar()
                title(lab[e])

                if ~isempty(ang_mf2)
                    # Grafique lineas que indican las direcciones principales 
                    norma = 2
            
                    for ang2 in ang_mf2
                        quiver(xnod[:,X],xnod[:,Y],                  # En el nodo grafique una línea
                                norma.*cos.(ang2), norma.*sin.(ang2), # indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                end

            elseif e == 3
                subplot(233)
                tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "jet",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                colorbar()
                title(lab[e])

                if ~isempty(ang_mt)

                    # Grafique lineas que indican las direcciones principales 
                    norma = 2
            
                    for ang3 in ang_mt
                        quiver(xnod[:,X],xnod[:,Y],                 # En el nodo grafique una línea
                                norma.*cos.(ang3), norma.*sin.(ang3),# indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                           plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                end


            else 
                
                
            end
        end

    return
 
end
