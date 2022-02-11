# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/losas/Kirchhoff_Love

# Traduciendo a JULIA por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


#Función para dibujar los diferentes gráficos del EF MZC,
#para la losa de Kirchhoff_Love
#Tener instalada la librería de matplotlib (PYTHON-pip install matplotlib)
#además Pkg.add("PyPlot") en consola de JULIA

using PyPlot,PyCall

function plot_mom_Q_ang(xnod, mom_Q, ang_, lab)

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
       
    for i = 1:length(mom_Q)
        norma = 1
   
        im = plt.tripcolor(xnod[:, X] , xnod[:, Y], triangles, mom_Q[i],  cmap = "bwr",
                       shading = "gouraud")
   
        colorbar(im,shrink=0.72)
        xlim(0, 2); ylim(0, 4); tight_layout()
        plt.gca().set_aspect("equal", adjustable="box")
        title(lab[i])
   
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
        for e = 1:nef
           nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                  plt.plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
        end
       end
       return
end