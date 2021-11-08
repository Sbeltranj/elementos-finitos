using DifferentialEquations
using Plots


function dibujar_deformada(E,A,I,L,x1,x2,y1,y2,qxloc,qyloc,ae)

#algunas constantes para facilitar la interpretación del código:

X  = 1; Y  = 2; 
X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6; 
v_ = 1; t_ = 2; M_ = 3; V_ = 4; u_ = 5; fax_ = 6;
u1 = 1; v1 = 2; t1 = 3; u2 = 4; v2 = 5; t2   = 6;


function diffeq1!(dydx,y,x,p)
    #=
    Esta función dibuja el elemento de pórtico deformado junto con sus 
    respectivos diagramas de fuerza axial, fuerza cortante y momento flector.
    El diagrama de momento flector se grafica en el lado opuesto de la fibra
    a tracción
    PARAMETROS DE ENTRADA (junto con algunos ejemplos):
    A = area
    E = E
    I = Ix local
    (x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
    qxloc = lambda x : x**2 # carga en la dir. del eje x local (function handle)
    qyloc = lambda x : 0    # carga en la dir. del eje y local (function handle)
    qe = np.array(
         [ 0.01,        # U1, V1, M1 reacciones del nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,        # U2, V2, M2 reacciones del nodo 2 en coord. locales
           0.02,
          -0.07 ])
    ae = np.array(
         [ 0.01,        # u1, v1, t1 desplazamientos nodo 1 en coord. locales
          -0.01,
           0.04,
          -0.01,        # u2, v2, t2 desplazamientos nodo 2 en coord. locales
           0.02,
          -0.07 ])
    esc_def    = 10     # escalamiento de la deformada
    esc_faxial = 10     # escalamiento del diagrama de axiales
    esc_V      = 10     # escalamiento del diagrama de cortantes
    esc_M      = 10     # escalamiento del diagrama de momentos 
    '''
    =#

    # aquí se implementa la ecuación diferencial para vigas de material
    # homogéneo y sección transversal constante (A, E, I, qxloc, qyloc las 
    # provee la función exterior)
    #      d^4 v(x)
    # E I ---------- = qyloc(x)
    #        dx^4
    #
    #      d^2 u(x)
    # A E ---------- = -qxloc(x)
    #        dx^2

    #            y[v_,:]          = v
    dydx[v_]   = y[t_]          # = theta
    dydx[t_]   = y[M_]/(E*I)    # = M/(EI)
    dydx[M_]   = y[V_]          # = V
    dydx[V_]   = qyloc(x)       # = qyloc
    dydx[u_]   = y[fax_]/(A*E)  # = u
    dydx[fax_] = qxloc(x)       # = faxial

end


#%% se definen las condiciones de frontera de la ecuacion diferencial 
function bc2!(residual,y,x,p)
   
    #residual == res
    #se específica la condición de frontera 
    #[1]   == apoyo izq   o condición inicial YL
    #[end] == apoyo derecho o condición final YR
    
    residual[1] = y[1][u_]     - ae[u1]  # y[1]   = YL
    residual[2] = y[1][v_]     - ae[v1]  # y[1]   = YL
    residual[3] = y[1][t_]     - ae[t1]  # y[1]   = YL
    residual[4] = y[end][u_]   - ae[u2]  # y[end] = YR
    residual[5] = y[end][v_]   - ae[v2]  # y[end] = YR
    residual[6] = y[end][t_]   - ae[t2]  # y[end] = YR

end

#BVProblem(function,conditions,solución tentativa EDO,tiempo-espacio sol EDO)
bvp = BVProblem(diffeq1!, bc2!, [0,0,0,0,0,0], (0.0,L))

#solve(bvp,algoritmo de solución, delta)
sol = solve(bvp, GeneralMIRK4(), dt=0.05)

#convertimos el vector tupla de la solución, en un array de 5x6 con el comando hcat:
sol1 = hcat(sol.u...)'

#%% Calculos intermedios
s     = sol.t

axial    = sol1[:,6]        # Fuerza axial [kN]
cortante = sol1[:,4]        # Fuerza cortante [kN]
momento  = sol1[:,3]        # Momento flector [kN/m]
u        = sol1[:,5]        # Desplazamiento horizontal de la viga [m]
v        = sol1[:,1]        # Desplazamiento vertical de la viga [m]

#%% rotacion de la solucion antes de dibujar
ang = atan(y2-y1, x2-x1)

T   = [ cos(ang)  -sin(ang)    # matriz de rotacion
        sin(ang)   cos(ang) ]


#%% Dibujar de deformada

pos = T*[ s +  10*u, 10*v ]; # escalamiento del diagrama


xx =  x1 .+ hcat(pos[1,:]...)
yy =  y1 .+ hcat(pos[2,:]...)

p1 = plot(1)
p1 = plot(xx,yy)
p1 = plot!([x1, x2], [y1, y2])


#%% Dibujar los diagramas de fuerza axial 
pos = T*[ s,   -10*axial ]; # escalamiento del diagrama

ss = x1 .+ hcat(pos[1,:]...)
aa = y1 .+ hcat(pos[2,:]...)



p2 = plot([x1, x2], [y1, y2])
p2 = plot!([x1; ss; x2], [y1; aa; y2])


#%% Dibujar los diagramas de fuerza cortante
pos = T*[ s,  cortante]; # escalamiento del diagrama

ss = x1 .+ hcat(pos[1,:]...)
vv = y1 .+ hcat(pos[2,:]...)


p3 = plot([x1, x2], [y1, y2])
p3 = plot!([x1; ss; x2], [y1; vv; y2])



#%% Dibujar los diagramas de momento flector
pos = T*[ s, momento ]; # escalamiento del diagrama

ss = x1 .+ hcat(pos[1,:]...)
mm = y1 .+ hcat(pos[2,:]...)


p4 = plot([x1, x2], [y1, y2])
p4 = plot!([x1; ss; x2], [y1; mm; y2])



end
