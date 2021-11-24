#Programa elaborado en JULIA 1.6.3
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
# -----------------------------------------------------------------
# Se usaron tres elementos isoparametricos lagrangianos cuadráticos
# -----------------------------------------------------------------
#Cargamos paquetes:

using Plots
using Polynomials


## defino las variables

E = 200e9;    # Pa               # módulo de elasticidad de la barra
A = (0.01)^2; # m^2              # área transversal de la barra
L = 2                            # m  # longitud de la barra
bb(x) = x.^2 - 2*x   # N/m       # fuerza axial aplicada sobre cada EF
P = 250                          # carga nodal al final de la barra
 
nef  = 1                         # nÚmero de elementos finitos (EF)
nno  = 3*nef + 1;                # nÚmero de nodos
ngdl = nno;
xnod = LinRange(0,L,nno)         # posición de los nodos

le = repeat([L/nef],nef)         # longitud de cada EF

LaG = [1 2 3 4];                 # definición de EFs con respecto a nodos

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

n_int_gl = 3;                 # orden de la cuadratura de Gauss-Legendre


x_gl, w_gl = gausslegendre_quad(n_int_gl)

# calcula las raices (xi_gl) y los pesos (w_gl) de polinomios de Legendre

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

## matriz de funciones de forma
Nforma(xi) = [ -(9*xi.^3)/16 + (9*xi.^2)/16 + xi/16 - 1/16    
                (27*xi.^3)/16 - (9*xi.^2)/16 - (27*xi)/16 + 9/16
                (27*xi)/16 - (9*xi.^2)/16 - (27*xi.^3)/16 + 9/16
                (9*xi.^3)/16 + (9*xi.^2)/16 - xi/16 - 1/16  ]

Bmat(Je, xi) = (1/Je).*[ (9*xi)/8 - (27*xi.^2)/16 + 1/16   
                (81*xi.^2)/16 - (9*xi)/8 - 27/16  
                27/16 - (81*xi.^2)/16 - (9*xi)/8  
                (27*xi.^2)/16 + (9*xi)/8 - 1/16  ];

### ensamblo la matriz de rigidez global y el vector de fuerzas nodales
##  equivalentes global
K  = zeros(ngdl,ngdl);   # matriz de rigidez global
De = E*A                # matriz constitutiva del elemento 
xe = Vector{Vector{Float64}}(undef, nef) # interpolación de la geometría

for e = 1:nef      # ciclo sobre todos los elementos finitos

    local idx, Je, Ke, fe
    idx  = LaG[e,:]
    
    
    Je = le[e]/2   # Jacobiano del elemento ( = dx_dxi)
    # vector de posiciones de nodos locales # = xnod[LaG[e,:]]';
    xe[e] = [xnod[LaG[e,1]], xnod[LaG[e,2]], xnod[LaG[e,3]], xnod[LaG[e,4]]];  
    
    # Cálculo de las matrices de rigidez y el vector de fuerzas nodales 
    # equivalentes del elemento
    Ke = zeros(4,4)
    fe = zeros(4,1)

    for m = 1:n_int_gl
       # matriz de deformacion del elemento
       local xi, Be
       xi = x_gl[m];
       Be = Bmat(Je, xi)'
       Ke += w_gl[m].*Be'*De.*Be.*Je; # matriz de rigidez del elemento e
 
       # vector de fuerzas nodales equivalentes
       local N, x_xi
       N = Nforma(xi)  # matr. de func. de forma
       x_xi = N'*xe[e] # interpola sobre la geometría (coord naturales a geométricas)
       fe .+=  (w_gl[m]*N'*bb(x_xi)*Je)'; # vector de fuerzas nodales equivalentes
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
nint = 200              # numero de puntos donde se interpolará dentro del EF
xi = collect(LinRange(-1,1,nint)) # coordenadas naturales

# matriz de funciones de forma
N = Nforma.(xi)
N = hcat(N...)'
xx    = Vector{Any}(undef,nef) # interpol de posiciones (geometria) en el elemento
uu    = Vector{Any}(undef,nef) # interpol desplazamientos en el elemento
axial = Vector{Any}(undef,nef) # fuerzas axiales en el elemento

for e in 1:nef       # ciclo sobre todas los elementos finitos

    local Je, Be, ae

    Je = le[e]/2      # Jacobiano del elemento ( = dx_dxi)
    Be = Bmat.(Je, xi) # matriz de deformacion del elemento
    Be = hcat(Be...)'
    # vector de desplazamientos nodales del elemento a^{(e)}
    ae = [ a[LaG[e,1]]
           a[LaG[e,2]]
           a[LaG[e,3]] 
           a[LaG[e,4]]] # = a(LaG(e,:))';
 
    # vector de posiciones de nodos locales
    
    xx[e] = N*xe[e] # interpola sobre la geometría (coord naturales a geométricas)
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
uexacto(x) = (x.*(4*L^3 .- 12*L^2 .- x.^3 + 4*x.^2 .+ 12*P))/(12*E*A) # solución análitica
x = collect(LinRange(0,L,200))               # 200 puntos unif/ distrib. entre 0 y L


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
Nexacta(x) =  E*A*((4*L^3 .- 12*L^2 .- x.^3 + 4*x.^2 .+ 12*P)/(12*E*A) + (x.*(-3*x.^2 + 8*x))/(12*E*A))       # solución análitica para la carga axial

fig_axial = plot(size = (80, 40))
fig_axial = scatter(x, Nexacta, label = "Solución analítica",legend=:bottomright,
                    title = "Comparación de la solución analÍtica con el MEF para la carga axial ",
                    titlefont=font(10,"Computer Modern"),
                    color = :blue,   
                    xaxis = "Eje X (m)",                         #nombre al eje x
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
savefig("./test.pdf")

display(fig_despla)
display(fig_axial)

#Fin