# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/1D/EF_barra_2_nodos

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

using Plots

#Package plotlyjs, (Pkg.add("plotlyjs")) en consola de Julia

#%% definición del problema
#% Calcule los desplazamientos y las reacciones en el empotramiento
#% de la barra mostrada
#%
#% | b (carga distribuida de magnitud b)
#% |->->->->->->->->->->->->->->->->
#% |====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
#% 1    2    3    4          nno-1  nno
#% |<----longitud L de la barra---->|   el área transversal de la barra es A

##  defino las variables
nef  = 3                      # número de elementos finitos (EF)
nno  = nef+1                  # número de nodos
ngdl = nno                    # número de grados de libertad
E    = 200e9    # Pa          # módulo de elasticidad de la barra
A    = (0.01)^2 # m^2         # área transversal de la barra
L    = 2        # m           # longitud de la barra
b    = 1000     # N/m         # fuerza axial aplicada sobre cada EF
P    = 250      # N           # carga nodal al final de la barra

xnod = LinRange(0, 2, 4)    # posición de los nodos
Le   = diff(xnod)           # longitud de cada EF (= repmat(L/nef, nef, 1))
k    = E*A./Le              # rigidez de cada EF

LaG = [(1:(nno-1)) (2:nno)] # definición de EFs con respecto a nodos

#%% Relación de cargas puntuales

f = zeros(ngdl,1); # vector de fuerzas nodales equivalentes global
f[nno] += P;        # relaciono la carga puntual en el nodo "nno"

#%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#%  equivalentes global
K = zeros(ngdl,ngdl)   # matriz de rigidez global

for e = 1:nef # ciclo sobre todos los elementos finitos
   idx         = LaG[e,:]
   K[idx,idx] += k[e]*[1 -1; -1 1];
   f[idx,:]   += ((b*Le[e])/2)*[1; 1];
end

#%% grados de libertad del desplazamiento conocidos y desconocidos
c = [1];    d = setdiff(1:ngdl, c);

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |
#|    | = |         ||    | - |    |     Recuerde que qc=0 (siempre)
#| qc |   | Kdc Kdd || ad |   | fc |

#%% extraigo las submatrices y específico las cantidades conocidas
Kcc = K[c,c]; Kdc = K[d,c]; fc = f[d]
Kdd = K[d,d]; Kcd = K[c,d]; fd = f[c];

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = 0;               # desplazamientos conocidos (en el gdl 1)

#%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac)
qd = Kcc*ac + Kcd*ad -fd;

a = zeros(ngdl,1);  q = zeros(ngdl,1);  # separo la memoria
a[c] .= ac;       q[c] = qd;
a[d] = ad;       #q[d] = qc = 0

#%% cálculo las cargas axiales en cada elemento finito

faxial = Array{Array{Float64}}(undef, nef,1) #separó memoria

for e = 1:nef # ciclo sobre todas los elementos finitos
   Be = [-1/Le[e] 1/Le[e]];
   ae = a[LaG[e,:]];
   faxial[e] = (E*A)*Be*ae; # = D*B(e)*a(e)
end
#%% imprimo los resultados

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Desplazamientos (m) = ")
display(a)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Fuerzas nodales equivalentes(N) =")
display(f)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Fuerzas nodales de equilibrio (N) = ")
display(q)

println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("Cargas axiales en cada elemento finito (N) =  ")
display(faxial)

#%% Gráfico la solución analítica y la solucion por el MEF
#% 1) gráfico los desplazamientos de la barra
u(x) = (-b*x.^2/2 + (P + b*L)*x)/(E*A); # solucion analítica para el despl.

xx = LinRange(0,L,100);                 # 100 puntos equidistantes entre 0 y L

#%% Nombró variable fig_despla  para ambos gráficos en un lienzo.

gr() #Pkg.add("GR")
#plotlyjs()  #instalar paquete plotlyjs (mejor interacción gráfica):
             #  Pkg.add("plotlyjs")
             #alternativamente comentar plotlyjs
# gráfico solución analítica
#plot! == hold on MatLab

fig_despla = plot(xx, u, label = "Solución analítica",
                    title = "Comparación de la solución analÍtica con el MEF para el desplazamiento ",
                    titlefont=font(10,"Computer Modern"),
                    color = :blue)

# gráfico solución por FEM
fig_despla  = plot!(xnod, a, label = "Solución por el FEM",
                    xaxis = "Eje X (m)",            #nombre al eje x
                    yaxis = "Desplazamiento (m)",   #nombre al eje y
                    color = :red)

#%% 2) gráfico la carga axial de la barra
faxial_exacta(x) =  (P + b*(L-x)) # solución analítica para la carga axial

#separó memoria x_nod, f_axial para solución por el FEM

x_nod     = Array{Array{Float64}}(undef, nef,1) #eje x gráfica
f_axial   = Array{Array{Float64}}(undef, nef,1) #eje Y gráfica

global fig_axi

fig_axi =plot() # creo lienzo en blanco para el ciclo for


# dibujo faxial exacta
fig_axi =  plot!(xnod, faxial_exacta,
           label = "Solución exacta",
           color = :red)

for e = 1:nef # ciclo sobre todas los elementos finitos

        x_nod[e]   = [xnod[e] xnod[e+1]]'     #eje x gráfica (xnod)
        f_axial[e] = [faxial[e] faxial[e]]'   #eje Y gráfica (faxial)

        # dibujo faxial por FEM
        if e  == 1
           fig_axi = plot!(x_nod[e], f_axial[e],
                     label = "Solución por el FEM", #etiqueta figura
                     xaxis = "Eje X (m)",           #nombre al eje x
                     yaxis = "Carga axial (N)",     #nombre al eje y
                     title = "Comparación de la solución analítica con el MEF para la fuerza axial ",
                     titlefont=font(10,"Computer Modern"),
                     color = :blue, shape = :circle)
         else
            fig_axi = plot!(x_nod[e], f_axial[e],
                      label = nothing,               #etiqueta figura
                      xaxis = "Eje X (m)",           #nombre al eje x
                      yaxis = "Carga axial (N)",     #nombre al eje y
                      title = "Comparación de la solución analítica con el MEF para la fuerza axial ",
                      titlefont=font(10,"Computer Modern"),
                      color = :blue, shape = :circle)
         end

end

#Imprimo las dos graficos, en un subplot (layout (2,1))
plot(fig_despla, fig_axi, layout=(2,1))

#%% Fin
