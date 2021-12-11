#= function  WoodArmer(Mx, My, Mxy)
    # Calculo los momentos de flexion Mxast_sup y Myast_sup 
    # asociados al refuerzo en la parte superior de la losa:
    # Se aplican las ecuaciones /*@\eqref{eq:Mxast_Myast_caso1_WA}@*/
    Mxast_sup = Mx + abs.(Mxy);
    Myast_sup = My + abs.(Mxy);
    
    if all(Mxast_sup .< 0) && all(Myast_sup .< 0)
        # no se requiere refuerzo en la parte superior de la losa
        Mxast_sup = 0;
        Myast_sup = 0;
    else
        ifelse.(Mxast_sup.<0, Mxast_sup, 0.)
            #Mxast_sup = 0;         #if all(Mxast_sup .< 0)
            #Myast_sup = My + abs.(Mxy.^2 ./Mx);
            if all(Myast_sup .< 0)
                Myast_sup = 0;
            end
           
        ifelse.(Myast_sup.<0, Myast_sup, 0.)
        #ifelse.(Myast_sup.<0, Mxast_sup, Mx + abs.(Mxy.^2 ./My))
        ifelse.(Mxast_sup.<0, Mxast_sup, 0.)
        
#=             #Mxast_sup = Mx + abs.(Mxy.^2 ./My);
            if all(Mxast_sup .< 0)
                Mxast_sup = 0;
            end     =#    
        
    end 
    
    # Calculo los momentos de flexion Mxast_inf y Myast_inf 
    # asociados al refuerzo en la parte inferior de la losa:
    # Se aplican las ecuaciones /*@\eqref{eq:Mxast_Myast_caso2_WA}@*/
    Mxast_inf = Mx - abs.(Mxy);
    Myast_inf = My - abs.(Mxy);
    if all(Mxast_inf .> 0) && all(Myast_inf .> 0)
        # no se requiere refuerzo en la parte inferior de la losa
        Mxast_inf = 0;
        Myast_inf = 0;
    else
        if all(Mxast_inf .> 0)
             Mxast_inf = 0;
             Myast_inf = My - abs.(Mxy.^2/Mx);
            if Myast_inf .> 0
                Myast_inf = 0;
            end        
        end    
        if all(Myast_inf .> 0)
            Mxast_inf = Mx - abs.(Mxy.^2/My);
            Myast_inf = 0;
            if all(Mxast_inf .> 0)
                Mxast_inf = 0;
            end                
        end
    end
    
    return Mxast_sup, Myast_sup, Mxast_inf, Myast_inf;
end
function dibujar_wood_armer(xnod, wood_A, lab)
    X,Y = 1,2
    
    NL1, NL2, NL3, NL4 = 1,2,3,4
    triangles = Vector{Vector{Int64}}(undef, 2*nef)
    for e = 1:nef
        # se arma la matriz de correspondencia (LaG) de la malla
        triangles[2*e - 1] = LaG[e, [NL1, NL2, NL4]] .- 1
        triangles[2*e - 0] = LaG[e, [NL2, NL3, NL4]] .- 1
    end
    fig, ax = subplots()
    for e = 1:length(wood_A)
            
        if e == 1
            subplot(141)
            tripcolor(xnod[:, X] , xnod[:, Y], triangles, wood_A[e],  cmap = "jet",
                      shading = "gouraud")
            xlim(0, 2); ylim(0, 4); tight_layout()
            plt.gca().set_aspect("equal", adjustable="box")
            colorbar()
            title(lab[e])
            for e = 1:nef
                # se dibujan las aristas
                nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                       plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
            end
        elseif e == 2
            subplot(142)
            tripcolor(xnod[:, X] , xnod[:, Y], triangles, wood_A[e],  cmap = "jet",
                      shading = "gouraud")
            xlim(0, 2); ylim(0, 4); tight_layout()
            plt.gca().set_aspect("equal", adjustable="box")
            colorbar()
            title(lab[e])
            for e = 1:nef
                # se dibujan las aristas
                nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                       plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
            end
        elseif e == 3
            subplot(143)
            tripcolor(xnod[:, X] , xnod[:, Y], triangles, wood_A[e],  cmap = "jet",
                      shading = "gouraud")
            xlim(0, 2); ylim(0, 4); tight_layout()
            plt.gca().set_aspect("equal", adjustable="box")
            colorbar()
            title(lab[e])
            for e = 1:nef
                # se dibujan las aristas
                nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                       plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
            end
        elseif e == 4
            subplot(144)
            tripcolor(xnod[:, X] , xnod[:, Y], triangles, wood_A[e],  cmap = "jet",
                      shading = "gouraud")
            xlim(0, 2); ylim(0, 4); tight_layout()
            plt.gca().set_aspect("equal", adjustable="box")
            colorbar()
            title(lab[e])
            for e = 1:nef
                # se dibujan las aristas
                nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL1]]
                       plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.15, color = "gray")
            end
        else
        end
    end
#=     ax.set_aspect("equal")
    ax.autoscale(tight=true)
    tight_layout() =#
    return
end =#