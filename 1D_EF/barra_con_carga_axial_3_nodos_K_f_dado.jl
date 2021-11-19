#Julia 1.6.3

#JULIA 1.6.3
## DEFINICIÓN DEL PROBLEMA
#=
Calcule los desplazamientos y las reacciones en el empotramiento
de la viga mostrada
| b (carga distribuida de magnitud b)
|->->->->->->->->->->->->->->->->
|====*====*====*====....====*====o-> P (carga puntual P en nodo nno)
1    2    3    4          nno-1  nno
|<----longitud L de la barra---->|   el área transversal de la barra es A
=#


# PROGRAMA ELABORADO POR: 
# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co

# Actualizando la versión 0.5.1 a 1.6.3
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

#Cargamos paquetes:

using Plots


# -----------------------------------------------------------------
# Se usaron tres elementos isoparametricos lagrangianos cuadraticos
# -----------------------------------------------------------------

## defino las variables

E = 200e9;    # Pa               # módulo de elasticidad de la barra
A = (0.01)^2; # m^2              # área transversal de la barra
L = 2                            # m                # longitud de la barra
b = 1000                         # N/m              # fuerza axial aplicada sobre cada EF
P = 250                          # carga nodal al final de la barra
 
nef  = 3                         # nÚmero de elementos finitos (EF)
nno  = 2*nef + 1;                # nÚmero de nodos
ngdl = nno;
xnod = LinRange(0,L,nno)         # posición de los nodos

le = repeat([L/nef],nef)         # longitud de cada EF

LaG = [1 2 3                     # definición de EFs con respecto a nodos
       3 4 5
       5 6 7];

## Relacion de cargas puntuales
f      = zeros(ngdl)  # vector de fuerzas nodales equivalentes global
f[nno] = P            # relaciono la carga puntual en el nodo "nno"

#%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#%  equivalentes global
K = zeros(ngdl,ngdl);   # matriz de rigidez global
De = E*A                # matriz constitutiva del elemento 

for e = 1:nef      # ciclo sobre todos los elementos finitos

    local idx
    idx  = LaG[e,:]
 
    Ke = (A*E/(3*le[e]))*[ 7    -8    1   # matriz de rigidez
                          -8    16   -8   # del elemento e
                           1    -8    7];

    fe = (b*le[e]/6)*[1; 4; 1]; # vector de fuerzas nodales equivalentes
 
    K[idx,idx] += Ke
    f[idx]     += fe
 end

## grados de libertad del desplazamiento conocidos y desconocidos
c = [ 1 ];    d = collect(2:nno)

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |  # recuerde que qc=0 (siempre)
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

## extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,c]; Kcd = K[c,d]; fd = f[c]
Kdc = K[d,c]; Kdd = K[d,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = [ 0 ];             # desplazamientos conocidos

## resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac)       # calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd    # calculo fuerzas de equilibrio desconocidas
a = zeros(nno);  a[c] = ac;  a[d] = ad # desplazamientos
q = zeros(nno);  q[c] = qd             # fuerzas nodales equivalentes


## se realizan unos calculos intermedios que necesitaremos mas adelante
nint = 10                # numero de puntos donde se interpolará dentro del EF
xi = collect(LinRange(-1,1,nint)) # coordenadas naturales

# matriz de funciones de forma
N = [xi.*(xi.-1)/2   (1 .+xi).*(1 .-xi)   xi.*(xi.+1)/2]
xx    = Vector{Any}(undef,nef) # interpol de posiciones (geometria) en el elemento
uu    = Vector{Any}(undef,nef) # interpol desplazamientos en el elemento
axial = Vector{Any}(undef,nef) # fuerzas axiales en el elemento

for e in 1:nef       # ciclo sobre todas los elementos finitos

    Je = le[e]/2      # Jacobiano del elemento ( = dx_dxi)
    Be = (1/Je)*[xi.-1/2  -2*xi  xi.+1/2] # matriz de deformacion del elemento
 
    # vector de desplazamientos nodales del elemento a^{(e)}
    ae = [ a[LaG[e,1]]
           a[LaG[e,2]]
           a[LaG[e,3]] ] # = a(LaG(e,:))';
 
    # vector de posiciones de nodos locales
    xe = [ xnod[LaG[e,1]]
           xnod[LaG[e,2]]
           xnod[LaG[e,3]] ] # = xnod(LaG(e,:))'
    xx[e] = N*xe # interpola sobre la geometría (coord naturales a geométricas)
    uu[e] = N*ae # interpola sobre los desplazamientos
 
    axial[e] = De*Be*ae # fuerzas axiales en elemento finito e
 end

## imprimo los resultados
# format short g
println("Desplazamientos (m) = ", a)
println("Fuerzas nodales equivalentes (N) = ", f)
println("Fuerzas nodales de equilibrio (N) = ", q)

## Grafico la solución análitica y la solución por el MEF
## 1) grafico los desplazamientos de la barra
uexacto(x) = (-b*x.^2/2 + (P + b*L)*x)/(E*A) # solución análitica
x = collect(LinRange(0,L,30))             # 30 puntos unif/ distrib. entre 0 y L


fig_despla = plot()
fig_despla = scatter(x, uexacto, label = "Solución analítica",
                    title = "Comparación de la solución analÍtica con el MEF para el desplazamiento ",
                    titlefont=font(10,"Computer Modern"),
                    color = :blue, 
                    xaxis = "Eje X (m)",              #nombre al eje x
                    yaxis = "Desplazamiento (m)", shape = :star5)   #nombre al eje y)

# gráfico solución por FEM
for e = 1:nef   # ciclo sobre todos los EFs
    if e == 1
        local fig_despla
        fig_despla  = plot!(xx[e], uu[e], color = :red, label = "Solución MEF") 
    else
        local fig_despla
        fig_despla  = plot!(xx[e], uu[e], color = :red, label = nothing) 
    end
                                                    
end

## 2) grafico la carga axial de la barra
Nexacta(x) = (P + b*(L-x))         # solución análitica para la carga axial

fig_axial = plot()
fig_axial = scatter(x, Nexacta, label = "Solución analítica",
                    title = "Comparación de la solución analÍtica con el MEF para la carga axial ",
                    titlefont=font(10,"Computer Modern"),
                    color = :blue, 
                    xaxis = "Eje X (m)",              #nombre al eje x
                    yaxis = "Carga axial (N)", shape = :star5)   #nombre al eje y)

# gráfico solución por FEM
for e = 1:nef   # ciclo sobre todos los EFs
    if e == 1
        local fig_axial
        fig_axial  = plot!(xx[e], axial[e], color = :red, label = "Solución MEF") 
    else
        local fig_axial
        fig_axial  = plot!(xx[e], axial[e], color = :red, label = nothing) 
    end
                                                    
end


display(fig_despla)
display(fig_axial)

#Fin
