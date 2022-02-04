# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/1D/EF_barra_n_nodos

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


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
#Cargamos paquetes:

using Plots
using Polynomials

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

function gausslegendre_quad(m)
# Integration using a Gauss-Legendre quadrature

## Calculation of the Legendre polynomials using Bonnet's recursion:
#           P_n(x) = ((2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x))/n

# Remember that  JULIA does not make 0-based indexing of arrays

P = Vector{Polynomial{Float64}}(undef, m+1)
P[0 + 1] =     Polynomial([1.0])       # P_{0}(x) = 1
P[1 + 1] = x = Polynomial([0.0, 1.0])  # P_{1}(x) = x

for n = 2:(m-1 + 1)
    P[n + 1] = ((2*n - 1)*x*P[n-1 + 1] - (n-1)*P[n-2 + 1])/n
end

    ## Roots
    xi = sort(roots(Polynomial(P[m+1])));

    ## Weights
    s = derivative(Polynomial(P[m+1]));

    w = 2.0 ./ ((1 .- xi.^2).*(s.(xi)).^2);

return xi, w
end

n_int_gl = 2;                 # orden de la cuadratura de Gauss-Legendre


x_gl, w_gl = gausslegendre_quad(n_int_gl)

# calcula las raíces (xi_gl) y los pesos (w_gl) de polinomios de Legendre

# >> [x_gl,w_gl] = gausslegendre_quad(1)
# x_gl = 0;
# w_gl = 2;
# >> [x_gl,w_gl] = gausslegendre_quad(2)
# x_gl = [  -0.577350269189626;  0.577350269189626 ];
# w_gl = [   1.000000000000000;  1.000000000000000 ];
# >> [x_gl,w_gl] = gausslegendre_quad(3)
# x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
# w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
# >> [x_gl,w_gl] = gausslegendre_quad(4)
# x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
# w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547; 0.347854845137453 ];

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
 
    Je = le[e]/2   # Jacobiano del elemento ( = dx_dxi)
    
    # Cálculo de las matrices de rigidez y el vector de fuerzas nodales 
    # equivalentes del elemento
    Ke = zeros(3,3)
    fe = zeros(3,1)

    for m = 1:n_int_gl
       # matriz de deformación del elemento
       local xi
       xi = x_gl[m];
       Be = (1/Je)*[xi-1/2 -2*xi xi+1/2];
       Ke += w_gl[m]*Be'*De*Be*Je; # matriz de rigidez del elemento e
 
       # vector de fuerzas nodales equivalentes
       local N
       N = [xi.*(xi-1)/2 (1+xi).*(1-xi) xi.*(xi+1)/2] # matr. de func. de forma
       fe .+=  w_gl[m]*N'*b*Je; # vector de fuerzas nodales equivalentes
    end
 
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
xx    = Vector{Any}(undef,nef) # interpol de posiciones (geometría) en el elemento
uu    = Vector{Any}(undef,nef) # interpol desplazamientos en el elemento
axial = Vector{Any}(undef,nef) # fuerzas axiales en el elemento

for e in 1:nef       # ciclo sobre todas los elementos finitos

    Je = le[e]/2      # Jacobiano del elemento ( = dx_dxi)
    Be = (1/Je)*[xi.-1/2  -2*xi  xi.+1/2] # matriz de deformación del elemento
 
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

## Grafico la solución analítica y la solución por el MEF
## 1) grafico los desplazamientos de la barra
uexacto(x) = (-b*x.^2/2 + (P + b*L)*x)/(E*A) # solución analítica
x = collect(LinRange(0,L,30))             # 30 puntos unif/ distrib. entre 0 y L


fig_despla = plot()
fig_despla = scatter(x, uexacto, label = "Solución analítica",
                    title = "Comparación de la solución analÍtica con el MEF para el desplazamiento ",
                    #titlefont=font(10,"Computer Modern"),
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
Nexacta(x) = (P + b*(L-x))         # solución analítica para la carga axial

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