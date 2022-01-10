# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/2D/modos_energia_nula

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co


using LinearAlgebra, PyPlot, PyCall, Printf
ENV["MPLBACKEND"]="qt5agg"
pygui(true)
close("all")          #cerrar ventanas

## Aquí se utiliza nullspace(K) para calcular el espacio nulo de la matrix.
## Programa para calcular los modos de energía nula de los EFs rectangulares
#  serendipitos de 4 y 8 nodos

X  = 1;
Y  = 2;
E  = 200;     # [GPa] módulo de elasticidad del elemento
nu = 0.33;    #       coeficiente de Poisson
t  = 0.10;    # [m]   espesor del elemento

## se selecciona el tipo de EF a analizar
nno = 4; # 4 u 8     # números de nodos del EF

if nno == 4
    ## Coordenadas del elemento
    xnod = [ 
    #  xi   eta     # nodo
       -1   -1      #  1
        1   -1      #  2
        1    1      #  3
       -1    1  ];  #  4

    ## Funciones de forma
    N(xi, eta) = [
                ((eta - 1)*(xi - 1))/4
                -((eta - 1)*(xi + 1))/4
                ((eta + 1)*(xi + 1))/4
                -((eta + 1)*(xi - 1))/4 ];

    ## Derivadas de N con respecto a xi
    dN_dxi(xi, eta) = [ 
                        eta/4 - 1/4              # dN1_dxi
                        1/4 - eta/4              # dN2_dxi
                        eta/4 + 1/4              # dN3_dxi
                       -(eta + 1)/4     ];      # dN4_dxi

    ## Derivadas de N con respecto a eta
    dN_deta(xi, eta) = [ 
                         xi/4 - 1/4               # dN1_deta
                        -(xi + 1)/4              # dN2_deta
                         xi/4 + 1/4               # dN3_deta
                         1/4 - xi/4       ];      # dN4_deta

    ## Modos de desplazamiento rigido
    a1 = [ 1 0 1 0 1 0 1 0 ]';
    a2 = [ 0 1 0 1 0 1 0 1 ]';

elseif nno == 8
    ## Coordenadas del elemento
    #  xi   eta     # nodo
    xnod = [ 
       -1   -1      #  1
        0   -1      #  2
        1   -1      #  3
        1    0      #  4
        1    1      #  5
        0    1      #  6
       -1    1      #  7
       -1    0  ];  #  8

    ## Funciones de forma N(xi,eta)
    N(xi, eta) = [ 
      -((eta - 1)*(xi - 1)*(eta + xi + 1))/4   # = N1
       ((xi^2 - 1)*(eta - 1))/2                # = N2
       ((eta - 1)*(xi + 1)*(eta - xi + 1))/4   # = N3
      -((eta^2 - 1)*(xi + 1))/2                # = N4
       ((eta + 1)*(xi + 1)*(eta + xi - 1))/4   # = N5
      -((xi^2 - 1)*(eta + 1))/2                # = N6
       ((eta + 1)*(xi - 1)*(xi - eta + 1))/4   # = N7
       ((eta^2 - 1)*(xi - 1))/2 ];             # = N8

    ## Derivadas de N con respecto a xi
    dN_dxi(xi, eta) =  [ 
       -((eta + 2*xi)*(eta - 1))/4              # dN1_dxi
       eta*xi - xi                              # dN2_dxi
        ((eta - 2*xi)*(eta - 1))/4              # dN3_dxi
       1/2 - eta^2/2                            # dN4_dxi
        ((eta + 2*xi)*(eta + 1))/4              # dN5_dxi
       -xi*(eta + 1)                            # dN6_dxi
       -((eta - 2*xi)*(eta + 1))/4              # dN7_dxi
       eta^2/2 - 1/2                        ];  # dN8_dxi

    ## Derivadas de N con respecto a eta
    dN_deta(xi, eta) = [ 
       -((2*eta + xi)*(xi - 1))/4               # dN1_deta
       xi^2/2 - 1/2                             # dN2_deta
       ((xi + 1)*(2*eta - xi))/4                # dN3_deta
       -eta*(xi + 1)                            # dN4_deta
       ((2*eta + xi)*(xi + 1))/4                # dN5_deta
       1/2 - xi^2/2                             # dN6_deta
       -((xi - 1)*(2*eta - xi))/4               # dN7_deta
       eta*(xi - 1)                         ];  # dN8_deta

    ## Modos de desplazamiento rigido
    a1 = [ 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 ]';   
    a2 = [ 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 ]';  
else 
    error("Solo se permite utilizar EFs de 4 y de 8 nodos")

end

## Se calculan y normalizan los modos de desplazamiento y rotación rígida
RmatT = [ cosd(45) -sind(45)
          sind(45)  cosd(45)]'; # matriz de rotación xr = Rmat*x


a1 = a1/norm(a1);
a2 = a2/norm(a2);
a3 = (xnod*RmatT)';   a3 = a3[:]/opnorm(a3)
mdrrr = [a1 a2 a3]; 

## Parámetros de la cuadratura de Gauss-Legendre
# se asumirá aquí el mismo orden de la cuadratura tanto en la dirección de
# xi como en la dirección de eta
n_gl = 2;   # orden de la cuadratura de Gauss-Legendre

include("gausslegendre_quad.jl")

x_gl, w_gl = gausslegendre_quad(n_gl);

## matriz constitutiva del elemento para TENSION PLANA
D = [ E/(1-nu^2)     E*nu/(1-nu^2)  0
      E*nu/(1-nu^2)  E/(1-nu^2)     0
      0              0              E/(2*(1+nu)) ];

# Calculo las matrices de rigidez y el vector de fuerzas nodales
# equivalentes del elemento
B = Array{Any}(undef,n_gl,n_gl) # contenedor para las matrices de deformación
K = zeros(2*nno, 2*nno);


for p = 1:n_gl
    for q = 1:n_gl

       xi_gl  = x_gl[p];
       eta_gl = x_gl[q];
       

       # Se evaluan las derivadas de las funciones de forma en los puntos
       # de integración de Gauss-Legendre
       xe = xnod[:,X];
       ye = xnod[:,Y];

       ddN_dxi  = dN_dxi(xi_gl, eta_gl)
       ddN_deta = dN_deta(xi_gl, eta_gl)
       
       dx_dxi  = sum(ddN_dxi .* xe);     dy_dxi  = sum(ddN_dxi  .* ye)
       dx_deta = sum(ddN_deta .* xe);    dy_deta = sum(ddN_deta .* ye)
       
       # Se ensambla la matriz Jacobiana del elemento
       Je = [ dx_dxi   dy_dxi
              dx_deta  dy_deta ]

       det_Je = zeros(n_gl,n_gl);
       det_Je[p,q] = det(Je)
       
       B[p,q] = zeros(3,2*nno);

      for i = 1:nno
         # Se ensambla la matriz de deformación del elemento B
         dNi_dx = (+dy_deta*ddN_dxi[i] - dy_dxi*ddN_deta[i])/det_Je[p,q]
         dNi_dy = (-dx_deta*ddN_dxi[i] + dx_dxi*ddN_deta[i])/det_Je[p,q]

         B[p,q][:,[2*i-1 2*i]] =         [dNi_dx      0      # aqui se ensambla
                                            0      dNi_dy    # y asigna la matriz
                                            dNi_dy   dNi_dx] # B_i
      end


       # se arma la matriz de rigidez del elemento e
       global K
       K += B[p,q]'*D*B[p,q]*det_Je[p,q]*t*w_gl[p]*w_gl[q];
    end
 end

## garantizar que K es simétrica, para evitar vectores propios complejos en
## eig
K = (K+K')/2; 

## Se calcula el espacio nulo de la matrix (contiene los vectores que generan dicho espacio)
null  = nullspace(K)
null_ = size(nullspace(K))

val, vec  = eigen(K)
#val = round.(val, digits = 3)
num_MEN = sum(val .< 1e-5)
idx = collect(1:1:null_[2])
Q, R = qr([mdrrr null[:,idx]]);
null[:,1:num_MEN] = Q[:,1:num_MEN]

## Se imprimen los vectores propios (recuerde que los modos de energia nula
## son aquellos para los cuales los valores propios son cero
##

xi  = [       collect(-1:0.1:1)'  +ones(1,length(-1:0.1:1))                collect(1:-0.1:-1)'   -ones(1,length(1:-0.1:-1)) ];
eta = [-ones(1,length(-1:0.1:1))               collect(-1:0.1:1)'    +ones(1,length(1:-0.1:-1))                collect(1:-0.1:-1)'   ];

modo = Array{Any}(undef,2*nno,1)

fig = plt.figure()
fig.suptitle("Puntos de integración = $(n_gl) x $(n_gl) MEN $(num_MEN)")

if     nno == 8 gs = fig.add_gridspec(ncols=2, nrows=2)
elseif nno == 4 gs = fig.add_gridspec(ncols=3, nrows=1)
end

for i = 1:null_[2]

    fig.add_subplot(gs[i])
    xlim(-2, 2); ylim(-2, 2);
    tight_layout()
    plt.gca().set_aspect("equal", adjustable="box")

    a = @sprintf("%.4E", val[i])
    title("λ$(i) = $(a)")
    modo[i] = xnod + reshape(null[:,i][1:nno*2],2,nno)';

    scatter(xnod[[ collect(1:nno)' 1],X],xnod[[ collect(1:nno)' 1],Y],
                label = nothing,color = :blue, marker = "*")
    plot(xi',eta',color = :blue, label = nothing, lw = 0.8)
    
    # Dibujar el elemento finito (utilizar la interpolación isoparamétrica)
    u = zeros(length(xi),1);
    v = zeros(length(xi),1);

    for j = 1:length(xi)
    u[j] = sum(N(xi[j],eta[j]).*modo[i][:,X]);
    v[j] = sum(N(xi[j],eta[j]).*modo[i][:,Y]);
    end

    plot(u,v,label = nothing, color="red", lw=0.8)

    scatter(modo[i][[ collect(1:nno)' 1],X]', modo[i][[ collect(1:nno)' 1],Y]',
            label = nothing, color = "red")

end