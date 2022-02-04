# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/losas

# Traduciendo a JULIA por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

using PyPlot,PyCall

function plot_mom_Q_ang(xnod, mom_Q, ang_, ang1, ang2, lab)

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
 
    fig = plt.figure()
    gs  = fig.add_gridspec(ncols=3, nrows=1)

    for i = 1:length(mom_Q)
     ax1 = fig.add_subplot(gs[1,i])
     im = ax1.tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[i],  cmap = "bwr",
                    shading = "gouraud")

     fig.colorbar(im,shrink=0.79)

     xlim(0, 2); ylim(0, 4); tight_layout()
     plt.gca().set_aspect("equal", adjustable="box")
     title(lab[i])

     for e = 1:nef
        nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
               plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
     end
      
     if i == 1
     if ~isempty(ang_)

        # Grafique lineas que indican las direcciones principales
        norma = 1
        for ang in ang_ #ang1 in ang

            plt.quiver(xnod[:,X],xnod[:,Y],                   # En el nodo grafique una línea
                    norma.*cos.(ang), norma.*sin.(ang), # indicando la dirección
                    headlength=0,
                    headwidth = 0,
                    headaxislength = 0,
                    scale = 8, pivot="middle")
        end
     end
    elseif i == 2
     if ~isempty(ang1)

        # Grafique lineas que indican las direcciones principales
        norma = 1
        for ang in ang1  #ang1 in ang
            plt.quiver(xnod[:,X],xnod[:,Y],                   # En el nodo grafique una línea
                    norma.*cos.(ang), norma.*sin.(ang), # indicando la dirección
                    headlength=0,
                    headwidth = 0,
                    headaxislength = 0,
                    scale = 8, pivot="middle")
        end
    end
    elseif i == 3
     if ~isempty(ang2)
        
        # Grafique lineas que indican las direcciones principales
        norma = 1
        for ang in ang2 #ang1 in ang
            plt.quiver(xnod[:,X],xnod[:,Y],                   # En el nodo grafique una línea
                    norma.*cos.(ang), norma.*sin.(ang), # indicando la dirección
                    headlength=0,
                    headwidth = 0,
                    headaxislength = 0,
                    scale = 8, pivot="middle")
        end 
    end
    end
    end
    # se dibujan las aristas

    

    return
end
