# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/losas

# Traduciendo a JULIA por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


#Función para dibujar los diferentes gráficos del EF QL9,

using PyPlot,PyCall

function plot_mom_Q_ang(xnod, mom_Q, ang_mf1, ang_mf2, ang_mt, lab)

    X,Y = 1,2
 
    NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8 = 1,2,3,4,5,6,7,8

    triangles = Vector{Vector{Int64}}(undef, 6*nef)


    for e = 1:nef

     # se arma la matriz de correspondencia (LaG) de la nueva malla triangular
        triangles[6*e - 5] = LaG[e, [NL1, NL2, NL8]] .- 1
        triangles[6*e - 4] = LaG[e, [NL2, NL3, NL4]] .- 1
        triangles[6*e - 3] = LaG[e, [NL4, NL5, NL6]] .- 1
        triangles[6*e - 2] = LaG[e, [NL2, NL4, NL6]] .- 1
        triangles[6*e - 1] = LaG[e, [NL2, NL6, NL8]] .- 1
        triangles[6*e - 0] = LaG[e, [NL6, NL7, NL8]] .- 1

    end
 
       
    fig, ax = subplots()

    for e = 1:length(mom_Q)
            
            if e == 1
                subplot(131)
                plt.tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "bwr",
                          shading = "gouraud")

                xlim(0, 2); ylim(0, 4); tight_layout()
                plt.gca().set_aspect("equal", adjustable="box")
                colorbar(shrink=0.79)
                title(lab[e])

                if ~isempty(ang_mf1)

                    # Grafique lineas que indican las direcciones principales
                    norma = 1

                    for ang1 in ang_mf1
                        plt.quiver(xnod[:,X],xnod[:,Y],                   # En el nodo grafique una línea
                                norma.*cos.(ang1), norma.*sin.(ang1), # indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end
               
                # se dibujan las aristas
                for e = 1:nef
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
                           plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                 end

            elseif e == 2
                subplot(132)
                
                plt.tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "bwr",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                plt.gca().set_aspect("equal", adjustable="box")
                colorbar(shrink=0.79)
                title(lab[e])

                if ~isempty(ang_mf2)
                    # Grafique lineas que indican las direcciones principales 
                    norma = 2
            
                    for ang2 in ang_mf2
                        plt.quiver(xnod[:,X],xnod[:,Y],                  # En el nodo grafique una línea
                                norma.*cos.(ang2), norma.*sin.(ang2), # indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
                    plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                end

            elseif e == 3
                subplot(133)
                plt.tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[e],  cmap = "bwr",
                          shading = "gouraud")
                xlim(0, 2); ylim(0, 4); tight_layout()
                plt.gca().set_aspect("equal", adjustable="box")
                colorbar(shrink=0.79)
                title(lab[e])

                if ~isempty(ang_mt)

                    # Grafique lineas que indican las direcciones principales 
                    norma = 2
            
                    for ang3 in ang_mt
                        plt.quiver(xnod[:,X],xnod[:,Y],                 # En el nodo grafique una línea
                                norma.*cos.(ang3), norma.*sin.(ang3),# indicando la dirección
                                headlength=0,
                                headwidth = 0,
                                headaxislength = 0,
                                scale = 8, pivot="middle")
                    end

                end

                for e = 1:nef
                    # se dibujan las aristas
                    nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
                    plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
                end
            else

            end

        end

    return
end