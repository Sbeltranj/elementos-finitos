#JULIA 1.7.1

# PROGRAMA ELABORADO POR: 
# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co

# Actualizando la versión 0.5.1 a 1.7.1
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

#Cargamos paquetes:

using Polynomials, PyPlot, LinearAlgebra, Statistics, SparseArrays, PyCall, WriteVTK

#Instale matplotlib de Python
ENV["MPLBACKEND"]="qt5agg"
pygui(true)
close("all")  # close all de MATLAB (función de PyPlot)

## ------------------------------------------------------------------------
## NOTA: este código SOLO es apropiado para tensión PLANA usando elementos
## rectangulares serendipitos de 8 nodos
## ------------------------------------------------------------------------

## definición del problema
# Calcule los desplazamientos y las reacciones en los empotramiento, las
# deformaciones y los esfuerzos dela estructura en tensión PLANA mostrada
# en la figura adjunta

include("gausslegendre_quad.jl")

function t2ft_R89(xnod, lado, carga, espesor)
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
 
    ## párametros de la cuadratura de Gauss-Legendre
    n_gl = 5                        # orden de la cuadratura de Gauss-Legendre
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
    te = zeros(16)
    te[[2*idx.-1; 2*idx]] = carga[:]
 
    ## Se calcula la integral
    suma   = zeros(16,16)
    N      = zeros(1,8)
    dN_dxi = zeros(1,8)

    for p = 1:n_gl
       N[idx] = NN(x_gl[p])
 
       matN = [ N[1] 0    N[2] 0    N[3] 0    N[4] 0    N[5] 0    N[6] 0    N[7] 0    N[8] 0
                0    N[1] 0    N[2] 0    N[3] 0    N[4] 0    N[5] 0    N[6] 0    N[7] 0    N[8] ]
 
       dN_dxi[idx] = dNN_dxi(x_gl[p])
 
       dx_dxi = dN_dxi*xnod[:,X]
       dy_dxi = dN_dxi*xnod[:,Y]
 
       # Toca poner ese [1] ya que dx_dxi y dy_dxi son matrices 1x1
       ds_dxi = hypot(dx_dxi[1], dy_dxi[1]) # sqrt(dx_dxi^2 + dy_dxi^2)
 
       suma += matN'*matN*ds_dxi*w_gl[p]
    end
 
    ft = espesor*suma*te
 
    return ft
 end


## defino las variables/constantes
X    = 1             # un par de constantes que ayudaran en la
Y    = 2             # lectura del código
Ee   = 200.0e9       # módulo de elasticidad del solido (Pa) = 200 GPa
nue  = 0.30          # coeficiente de Poisson
te   = 0.01          # espesor del solido (m)
rhoe = 7850.0        # densidad (kg/m^3)
g    = 9.81          # aceleración de la gravedad (m/s^2)
be = [0; -rhoe*g]    # vector de fuerzas másicas del elemento

MALLA = 3   # MALLA=1 gráfico, MALLA=2 la generada con ANSYS

## cargar
# xnod - posición de los nodos
# LaG  - definición de elementos finitos con respecto a nodos
if     MALLA == 1    include("malla1.jl")
elseif MALLA == 2    include("malla2.jl")
elseif MALLA == 3    include("malla3.jl")
else                 error("Malla no especificada")
end


nno  = size(xnod,1)  # número de nodos (número de filas de xnod)
ngdl = 2*nno         # número de grados de libertad (dos por nodo)
gdl  = [1:2:ngdl   2:2:ngdl] # nodos vs grados de libertad
nef = size(LaG,1)    # número de EFs (número de filas de LaG)

## Se definen las restricciones
ngdl_res = size(restricciones,1) # número de grados de libertad restringidos
restric  = zeros(ngdl_res, 2)#Array{Array{Float64}}(undef, ngdl_res,2)

for i = 1:ngdl_res
   #                    nodo                     dirección                desplazamiento
   restric[i,:] = [ gdl[Int(restricciones[i,1]), Int(restricciones[i,2])] restricciones[i,3] ]
end


f = zeros(ngdl)         # vector de fuerzas nodales equivalentes global
if MALLA == 1
   f[gdl[13,Y]] = -5000 # carga puntual en el nodo 13 dir Y
   f[gdl[21,Y]] = -5000 # carga puntual en el nodo 21 dir Y
elseif MALLA == 2
   f[gdl[23,Y]] = -5000 # carga puntual en el nodo 13 dir Y
   f[gdl[29,Y]] = -5000 # carga puntual en el nodo 21 dir Y
elseif MALLA == 3
   # No tenemos cargas puntuales en esta estructura
else
   error("Malla no especificada")
end

## Se definen las cargas distribuidas
if MALLA == 1
    carga_distr = Array{Float64}(undef,0,6)
elseif MALLA == 2
    carga_distr = Array{Float64}(undef,0,6)
elseif MALLA == 3
 #  elem lado  tix   tiy   tjx    tjy  tkx  tky
    carga_distr = [
    42   123   0         0 0     50000 0    100000  # 42 702 710 717 123
    43   123   0    100000 0    100000 0    100000  # 43 717 727 733 123
    44   123   0    100000 0    100000 0    100000  # 44 733 740 747 123
    45   123   0    100000 0    100000 0    100000  # 45 747 754 760 123
    46   123   0    100000 0    100000 0    100000  # 46 760 767 772 123
    47   123   0    100000 0    100000 0    100000  # 47 772 776 779 123
    48   123   0    100000 0    100000 0    100000  # 48 779 780 778 123
    49   123   0    100000 0    100000 0    100000  # 49 778 774 770 123
    50   123   0    100000 0    100000 0    100000  # 50 770 766 761 123
    51   123   0    100000 0     50000 0       0 ]  # 51 761 757 750 123
else
    error("Malla no especificada")
end

nlcd = size(carga_distr,1) # número de lados con carga distribuida



figure(1)
cg = zeros(nef, 2) # almacena el centro de gravedad
for e = 1:nef
   plt.plot(xnod[LaG[e,[1:8; 1]],X], xnod[LaG[e,[1:8; 1]],Y],
         color="k", linestyle="-")

   # Calculo la posición del centro de gravedad del triangulo
   cg[e,:] = [ mean(xnod[LaG[e,[1 3 5 7]],X]) mean(xnod[LaG[e,[1 3 5 7]],Y]) ]

   plt.text(cg[e, X], cg[e, Y], "$e", fontsize=5, color=[1,0,0],
         horizontalalignment="center", verticalalignment="center")
end
plot(xnod[:,X], xnod[:,Y], "b.")
# text(xnod[:,X], xnod[:,Y], num2str((1:nno)'), fontsize=16)
axis("equal") # falta tight
title("Malla de elementos finitos")

## Funciones de forma serendipitas del elemento rectangular de 8 nodos:
# NOTA estas funciones de forma y sus derivadas se encontraron con el
# programa c5_funciones_forma_lagrangianos_rect_2D_8_nodos.m (https://github.com/diegoandresalvarez)

Nforma(xi,eta) = [
      -((eta - 1)*(xi - 1)*(eta + xi + 1))/4       # N1
      ((xi^2 - 1)*(eta - 1))/2                     # N2
      ((eta - 1)*(xi + 1)*(eta - xi + 1))/4        # N3
      -((eta^2 - 1)*(xi + 1))/2                    # N4
      ((eta + 1)*(xi + 1)*(eta + xi - 1))/4        # N5
      -((xi^2 - 1)*(eta + 1))/2                    # N6
      ((eta + 1)*(xi - 1)*(xi - eta + 1))/4        # N7
      ((eta^2 - 1)*(xi - 1))/2                 ]   # N8

## Derivadas de N con respecto a xi
dN_dxi(xi,eta) = [
      -((eta + 2*xi)*(eta - 1))/4                  # dN1_dxi
      eta*xi - xi                                  # dN2_dxi
      ((eta - 2*xi)*(eta - 1))/4                   # dN3_dxi
      1/2 - eta^2/2                                # dN4_dxi
      ((eta + 2*xi)*(eta + 1))/4                   # dN5_dxi
      -xi*(eta + 1)                                # dN6_dxi
      -((eta - 2*xi)*(eta + 1))/4                  # dN7_dxi
      eta^2/2 - 1/2                            ]   # dN8_dxi

## Derivadas de N con respecto a eta
dN_deta(xi,eta) = [
      -((2*eta + xi)*(xi - 1))/4                   # dN1_deta
      xi^2/2 - 1/2                                 # dN2_deta
      ((xi + 1)*(2*eta - xi))/4                    # dN3_deta
      -eta*(xi + 1)                                # dN4_deta
      ((2*eta + xi)*(xi + 1))/4                    # dN5_deta
      1/2 - xi^2/2                                 # dN6_deta
      -((xi - 1)*(2*eta - xi))/4                   # dN7_deta
      eta*(xi - 1)                             ]   # dN8_deta

## parámetros de la cuadratura de Gauss-Legendre
# se asumirá aquí el mismo orden de la cuadratura tanto en la dirección de
# xi como en la dirección de eta
n_gl = 2    # orden de la cuadratura de Gauss-Legendre


# El comando:
x_gl, w_gl  = gausslegendre_quad(n_gl)

## ensamblo la matriz de rigidez global y el vector de fuerzas nodales
#  equivalentes global


K = spzeros(ngdl,ngdl)        # matriz de rigidez global como RALA (sparse)
N = Array{Any}(undef,nef,n_gl,n_gl) # contenedor para las matrices de forma
B = Array{Any}(undef,nef,n_gl,n_gl) # contenedor para las matrices de deformación

# matriz constitutiva del elemento para tensión PLANA
De = [ Ee/(1-nue^2)     Ee*nue/(1-nue^2)  0
       Ee*nue/(1-nue^2) Ee/(1-nue^2)      0
       0                0                 Ee/(2*(1+nue)) ]


idx   = Array{Array{Int64}}(undef, nef,1) 

for e = 1:nef  # ciclo sobre todos los elementos finitos
    idx[e] = [  gdl[LaG[e,1],:];  gdl[LaG[e,2],:];
                gdl[LaG[e,3],:];  gdl[LaG[e,4],:];
                gdl[LaG[e,5],:];  gdl[LaG[e,6],:];
                gdl[LaG[e,7],:];  gdl[LaG[e,8],:] ]

       # Calculo las matrices de rigidez y el vector de fuerzas nodales
   # equivalentes del elemento
   Ke = zeros(16,16)
   fe = zeros(16)
   det_Je = zeros(n_gl,n_gl) # en esta matriz se almacenaran los Jacobianos

   for p = 1:n_gl
      for q = 1:n_gl
         xi_gl  = x_gl[p]
         eta_gl = x_gl[q]

         # Se evalúan las funciones de forma en los puntos de integración
         # de Gauss-Legendre
         NNforma = Nforma(xi_gl, eta_gl)

         # Se evalúan las derivadas de las funciones de forma en los puntos
         # de integración de Gauss-Legendre
         ddN_dxi  = dN_dxi( xi_gl, eta_gl);     xe = xnod[LaG[e,:],X]
         ddN_deta = dN_deta(xi_gl, eta_gl);     ye = xnod[LaG[e,:],Y]

         dx_dxi  = sum(ddN_dxi  .* xe);      dy_dxi  = sum(ddN_dxi  .* ye)
         dx_deta = sum(ddN_deta .* xe);      dy_deta = sum(ddN_deta .* ye)

         # Se ensambla la matriz Jacobiana del elemento
         Je = [ dx_dxi   dy_dxi
                dx_deta  dy_deta ]

         # Se calcula el determinante del Jacobiano
         det_Je[p,q] = det(Je)

         N[e,p,q] = zeros(2,2*8)
         B[e,p,q] = zeros(3,2*8)

         for i = 1:8
            # Se ensambla la matriz de funciones de forma N
            N[e,p,q][:,[2*i-1 2*i]] = [ NNforma[i]  0
                                        0           NNforma[i] ]

            # Se ensambla la matriz de deformación del elemento B
            dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je[p,q]
            dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je[p,q]
            B[e,p,q][:,[2*i-1 2*i]] = [ dNi_dx 0          # aquí se ensambla
                                        0      dNi_dy     # y asigna la matriz
                                        dNi_dy dNi_dx ]   # B_i
         end

         # se arma la matriz de rigidez del elemento e
         Ke += B[e,p,q]'*De*B[e,p,q]*det_Je[p,q]*te*w_gl[p]*w_gl[q]

         # vector de fuerzas nodales equivalentes
         fe += N[e,p,q]'*be*det_Je[p,q]*te*w_gl[p]*w_gl[q]
      end
   end;
   if any(any(det_Je .<= 0))
    error("Existen elementos con det_Je negativo en el elemento $e.\n")
 end

 K[idx[e],idx[e]] += Ke
 f[idx[e],:]      += fe
end

## Relación de las cargas superficiales (vector ft)
ft = zeros(ngdl)   # fuerzas nodales equivalentes de cargas superficiales
for i = 1:nlcd
   e     = carga_distr[i,1]
   lado  = carga_distr[i,2]
   carga = carga_distr[i,3:8]
   fte   = t2ft_R89(xnod[LaG[e,1:8],[X, Y]], lado, carga, te)
   ft[idx[e]] += fte
end

# Agregó al vector de fuerzas nodales equivalentes las fuerzas
# superficiales calculadas
f += ft

## Muestro la configuración de la matriz K (K es rala)
figure(2)
spy(K)
title("Los puntos representan los elementos diferentes de cero")

## grados de libertad del desplazamiento conocidos y desconocidos
c = restric[:,1];   d = setdiff(1:ngdl,c)
c = round.(Int, c)  # se convierte en Int64
# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos

#| qd |   | Kcc Kcd || ac |   | fd |
#|    | = |         ||    | - |    |
#| qc |   | Kdc Kdd || ad |   | fc |

## extraigo las submatrices y especifico las cantidades conocidas
Kcc = K[c,c]; Kcd = K[c,d]; fd = f[c]
Kdc = K[d,c]; Kdd = K[d,d]; fc = f[d]

# f = vector de fuerzas nodales equivalentes
# q = vector de fuerzas nodales de equilibrio del elemento
# a = desplazamientos
ac = restric[:,2]    # desplazamientos conocidos
qc = zeros(size(d))  # cargas de equilibrio en nodos libres ( = 0 siempre)
qc = round.(Int, qc)
# es lo mismo que qc = zeros(size(d)); de MATLAB

## resuelvo el sistema de ecuaciones
ad = Kdd \ ((fc + qc)-Kdc*ac)     # calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd     # calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl);  a[c] = ac;  a[d] = ad   # desplazamientos
q = zeros(ngdl);  q[c] = qd;  q[d] = qc   # fuerzas nodales equivalentes

## imprimo los resultados
# format short g
## Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,2,nno)'

escala = 50000                # factor de escalamiento de la deformada
xdef = xnod + escala*delta    # posición de la deformada
figure(3)
for e = 1:nef
   plt.plot(xnod[LaG[e,[1:8;1]],X], xnod[LaG[e,[1:8;1]],Y], color="r", linestyle="-")   # original
   plt.plot(xdef[LaG[e,[1:8;1]],X], xdef[LaG[e,[1:8;1]],Y], color="b", linestyle="-")   # deformada
end
axis("tight")
axis("equal")
legend(["posición original","posición deformada"], loc="lower right")
title("Deformada escalada $escala veces")


## Se calcula para cada elemento las deformaciones y los esfuerzos
def = Array{Any}(undef,nef,n_gl,n_gl)
esf = Array{Any}(undef,nef,n_gl,n_gl)

for e = 1:nef
   ae = a[idx[e]]            # desplazamientos de los gdl del elemento e

   for pp = 1:n_gl
      for qq = 1:n_gl
         def[e,pp,qq] = B[e,pp,qq]*ae     # calculo las deformaciones
         esf[e,pp,qq] = De*def[e,pp,qq]   # calculo los esfuerzos
      end
   end
end

## Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno)  # número de elementos adyacentes
sx  = zeros(nno)
sy  = zeros(nno)
sz  = zeros(nno)
txy = zeros(nno)
txz = zeros(nno)
tyz = zeros(nno)

ex  = zeros(nno)
ey  = zeros(nno)
gxy = zeros(nno)

A = [
     3^(1/2)/2 + 1             -1/2             -1/2   1 - 3^(1/2)/2
   3^(1/2)/4 + 1/4  1/4 - 3^(1/2)/4  3^(1/2)/4 + 1/4 1/4 - 3^(1/2)/4
              -1/2    1 - 3^(1/2)/2    3^(1/2)/2 + 1            -1/2
   1/4 - 3^(1/2)/4  1/4 - 3^(1/2)/4  3^(1/2)/4 + 1/4 3^(1/2)/4 + 1/4
     1 - 3^(1/2)/2             -1/2             -1/2   3^(1/2)/2 + 1
   1/4 - 3^(1/2)/4  3^(1/2)/4 + 1/4  1/4 - 3^(1/2)/4 3^(1/2)/4 + 1/4
              -1/2    3^(1/2)/2 + 1    1 - 3^(1/2)/2            -1/2
   3^(1/2)/4 + 1/4  3^(1/2)/4 + 1/4  1/4 - 3^(1/2)/4 1/4 - 3^(1/2)/4 ];

for e = 1:nef
   sx[LaG[e,:],:] .+=  A * [ esf[e,1,1][1]
                                              esf[e,1,2][1]
                                              esf[e,2,1][1]
                                              esf[e,2,2][1] ]

   sy[LaG[e,:],:] .+=  A * [ esf[e,1,1][2]
                                              esf[e,1,2][2]
                                              esf[e,2,1][2]
                                              esf[e,2,2][2] ]

   txy[LaG[e,:],:] .+= A * [ esf[e,1,1][3]
                                              esf[e,1,2][3]
                                              esf[e,2,1][3]
                                              esf[e,2,2][3] ]

   ex[LaG[e,:],:] .+= A * [ def[e,1,1][1]
                                              def[e,1,2][1]
                                              def[e,2,1][1]
                                              def[e,2,2][1] ]

   ey[LaG[e,:],:] .+= A * [ def[e,1,1][2]
                                              def[e,1,2][2]
                                              def[e,2,1][2]
                                              def[e,2,2][2] ]

   gxy[LaG[e,:],:] .+=  A * [ def[e,1,1][3]
                                              def[e,1,2][3]
                                              def[e,2,1][3]
                                              def[e,2,2][3] ]

   num_elem_ady[LaG[e,:],:] .+= 1
end

## Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;  ex  =  ex./num_elem_ady
sy  =  sy./num_elem_ady;  ey  =  ey./num_elem_ady
txy = txy./num_elem_ady;  gxy = gxy./num_elem_ady

## Se calculan las deformación ez en tensión plana
ez  = -(nue/Ee)*(sx+sy)


# Ver: http://stackoverflow.com/questions/29443369/how-to-make-a-custom-colormap-using-pyplot-not-matplotlib-proper
@pyimport matplotlib.colors as my_colors

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

   #Actualizando a tripcolor:
   #https://matplotlib.org/stable/gallery/images_contours_and_fields/tripcolor_demo.html

      val_max = maximum(abs.(esfdef))
      fig, ax = subplots()
       # se grafica la malla de EFS, los colores en cada triángulo y las curvas 
       # de nivel
      for e = 1:nef
         # se dibujan las aristas
         nod_ef = LaG[e, [NL1, NL2, NL3, NL4, NL5, NL6, NL7, NL8, NL1]]
                plot(xnod[nod_ef, X], xnod[nod_ef, Y], lw = 0.5, color = "gray")
      end

      im = ax.tripcolor(xnod[:, X], xnod[:, Y], triangle, esfdef,  cmap = "jet",
                        shading = "gouraud", vmin = -val_max, vmax = val_max)

      ax.tricontour(xnod[:, X], xnod[:, Y], triangle, esfdef, 20)
      ax.tricontour(xnod[:, X], xnod[:, Y], triangle, esfdef, levels=[0], linewidths=3)

      fig.colorbar(im, ax = ax, format = "%6.3g")

      if ~isempty(angulos)
         # Grafique lineas que indican las direcciones principales de sigma_1
         norma = 1 # = esf si quiere proporcional
   
         for ang in angulos
            plt.quiver(xnod[:,X],xnod[:,Y],              # En el nodo grafique una línea
                   norma.*cos.(ang),norma.*sin.(ang),# indicando la dirección
                   headlength=0,
                   headwidth = 0,
                   headaxislength = 0,
                   pivot="middle")
         end
         # scatter(xnod[:,X],xnod[:,Y], color="k", s=1)  # para poner el punto
         # http://matplotlib.org/examples/pylab_examples/quiver_demo.html
      end
      #ylabel(lab)
      ax.set_xlabel("x [m]")
      ax.set_ylabel("y [m]")
      ax.set_title(lab, fontsize=20)
      axis("equal") # tight
   return

end

plot_def_esf_ang(xnod, ex,  [], L"\epsilon_x(x,y)")
plot_def_esf_ang(xnod, ey,  [], L"\epsilon_y(x,y)")
plot_def_esf_ang(xnod, ez,  [], L"\epsilon_z(x,y)")
plot_def_esf_ang(xnod, gxy, [], L"\gamma_{xy}(x,y)")

## Se imprimen y grafican los esfuerzos en los nodos
plot_def_esf_ang(xnod, sx,  [], L"\sigma_x(x,y) [Pa]")
plot_def_esf_ang(xnod, sy,  [], L"\sigma_y(x,y) [Pa]")
plot_def_esf_ang(xnod, txy, [], L"\tau_{xy}(x,y) [Pa]")

## Se calculan y grafican para cada elemento los esfuerzos principales y
## sus direcciones
# NOTA: esto solo es válido para el caso de TENSIÓN PLANA).
# En caso de DEFORMACIÓN PLANA se deben calcular los valores y vectores
# propios de la matriz de tensiónes de Cauchy
#   [dirppales{e}, esfppales{e}] = eig([sx  txy 0    # matriz de esfuerzos
#                                       txy sy  0    # de Cauchy
#                                       0   0   0]);
s1   = (sx+sy)/2 + sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal máximo
s2   = (sx+sy)/2 - sqrt.(((sx-sy)/2).^2+txy.^2) # esfuerzo normal mínimo
tmax = (s1-s2)/2                               # esfuerzo cortante máximo
ang  = 0.5*atan.(2*txy, sx-sy) # ángulo de inclinación de s1

## imprimo los resultados
## s1, s2, taumax
plot_def_esf_ang(xnod, s1,   [ang],                  L"\sigma_1(x,y) [Pa]")
plot_def_esf_ang(xnod, s2,   [ang.+pi/2],            L"\sigma_2(x,y) [Pa]")
plot_def_esf_ang(xnod, tmax, [ang.+pi/4, ang.-pi/4], L"\tau_{max}(x,y) [Pa]")

## Cálculo de los esfuerzos de von Mises
s3 = zeros(size(s1))   # s3 = zeros(size(s1)) de MATLAB
sv = sqrt.(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2)

#println("Nodo,Esfuerzos de von Mises (Pa) = ")
#println([(1:nno)'  sv]);
plot_def_esf_ang(xnod, sv, [], L"\sigma_v(x,y) [Pa]")
title("Esfuerzos de von Mises (Pa)")

#Consulte la documentación:
#https://jipolanco.github.io/WriteVTK.jl/dev/grids/unstructured/#Unstructured-grid

## se reportan resultados .vtu, para visualizar en  paraview

cells   = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,nef) 
for e = 1:nef

 ## Nota,de acuerdo al manual de Paraview el EF de 8 nodos se denomina VTK_QUADRATIC_QUAD (=23)
 ## Sin embargo al implementarlo en JULIA no es posible llevarlo a cabo, por esto en connectivity,
 ## He tomado los vertices del EF así como se muestra a continuación. ##Con HexaHedron (3D) parece funcionar...
 
 global cells[e] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON,  vec(LaG[e,[1 3 5 7 8 2 4 6]]) )
 
end

vtkfile = vtk_grid("Q8_element", xnod[:,X],xnod[:,Y], cells) 

vtkfile["uv"]  = a 

vtkfile["s_x"]   = sx;    vtkfile["s1"] = s1;      vtkfile["ex"] = ex; vtkfile["gxy"] = gxy;
vtkfile["s_y"]   = sy;    vtkfile["s2"] = s2;      vtkfile["ey"] = ey;
vtkfile["t_xy"]  = txy; vtkfile["Tmax"] = tmax;  vtkfile["ez"]   = ez;

vtkfile["sv"] = sv
vtkfile["n1"] = [cos.(ang)           sin.(ang)                    ]
vtkfile["n2"] = [cos.(ang .+ pi/2)           sin.(ang .+ pi/2)    ]

outfiles = vtk_save(vtkfile)

println("Se han reportado los datos en formato .vtu para ser visualizados en ParaView")