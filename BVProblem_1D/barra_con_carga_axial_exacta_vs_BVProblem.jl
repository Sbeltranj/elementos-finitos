#JULIA 1.6.3

#se cargan paquetes:
using DifferentialEquations
using Plots

## %% Definición del problema
# Calcule los desplazamientos y las reacciones en el empotramiento
# de la barra mostrada resolviendo la ecuación diferencial numéricamente con
# la funcion BVProblem
#
# | b (carga distribuida de magnitud b)
# |->->->->->->->->->->->->->->->->
# |================================o--> P (carga puntual P en extremo derecho)
# |<----longitud L de la barra---->|    el área transversal de la barra es A

E = 200e9     # Pa               # módulo de elasticidad de la barra
A = (0.01)^2  # m^2              # área transversal de la barra
L = 2         # m                # longitud de la barra
b = 1000      # N/m              # fuerza axial aplicada sobre cada EF
P = 250       # N                # carga nodal al final de la barra

#%% Solución de la ecuacion diferencial

# Solución numérica usando solve_bvp (boundary value problem)
#   d /           du(x)  \
# ----| E(x) A(x)------- | + b(x) en x \in [0,L]     dado u(0) = 0
#  dx \            dx    /                                faxial(L) = P


# En el caso mas general E, A y b son funciones. Escriba aquí las funciones
# como tal en caso de tener un caso mas general. 
EE(x) = E
AA(x) = A
bb(x) = b

#%% Se plantea el sistema de ecuaciones:
#consultar BVProblem: https://diffeq.sciml.ai/stable/tutorials/bvp_example/

#sistema de ecuaciones diferenciales acopladas:
function diffeq!(dy,y,p,x)
    dy[1] = y[2]/(EE(x)*AA(x))  # dy1/dx = y2/E(x)*A(x)
    dy[2] = -bb(x)              # dy2/dx = -b(x)
end

#Condición de contorno:
function bc!(residual,y,p,x)
   residual[1] = y[1][1] - 0    # y1(0) = 0
   residual[2] = y[end][2] - P  # y2(L) = P
end


#BVProblem(function,conditions,solución tentativa EDO,tspan)
bvp = BVProblem(diffeq!, bc!, [0,0], (0.0,L))

#solve(bvp,algoritmo de solución, delta)
sol = solve(bvp, GeneralMIRK4(), dt=0.05)

#%% Solución analítica

# 41 puntos uniformemente distrib. entre 0 y L
x = collect(LinRange(0,2,41))

u_exacta(x)      = (-b*x.^2/2 + (P + b*L)*x)/(E*A) #desplazamientos
faxial_exacta(x) = (P + b*(L-x))                   #fuerza axial           

#%% Se reportan los errores en el cálculo
error_en_u      = maximum(abs.(u_exacta(x)- first.(sol.u)))
println("Maximo error en el cálculo del desplazamiento =  $(error_en_u) m")

error_en_faxial = maximum(abs.(faxial_exacta.(x) - last.(sol.u)))
println("Maximo error en el cálculo de la fuerza axial =  $(error_en_faxial) N")


#%% Grafico la solución analítica y la solución por el la función BVProblem


# 1) grafico los desplazamientos de la barra

#se dibuja la función exacta en la variable fig_despla
fig_despla  =  plot(x, u_exacta, label = "Solución exacta", #etiqueta
                    xaxis = "Eje X (m)",                    #nombre al eje x
                    yaxis = "Desplazamiento (m)",           #nombre al eje y
                    color = :red)


#se dibuja la solución BVProblem en la variable fig_despla
fig_despla  =  scatter!(x, first.(sol.u),             # ! == hold on MATLAB 
                    label = "Solución por BVProblem", #etiqueta
                    xaxis = "Eje X (m)",              #nombre al eje x
                    yaxis = "Desplazamiento (m)",     #nombre al eje y
                    title = "Comparación de la solución analítica vs BVProblem para el desplazamiento",
                    titlefont=font(10,"Computer Modern"),
                    shape = :star5, color = :blue)


# 2) grafico la carga axial de la barra

fig_axial  =   plot(x, faxial_exacta, label = "Solución exacta", #etiqueta
                    xaxis = "Eje X (m)",                         #nombre al eje x
                    yaxis = "Desplazamiento (m)",                #nombre al eje y
                    color = :red)

#se dibuja la solución BVProblem
fig_axial  =  scatter!(x, last.(sol.u) , label = "Solución por BVProblem",
                        xaxis = "Eje X (m)",            #nombre al eje x
                        yaxis = "Carga axial (N)",      #nombre al eje y
                        title = "Comparación de la solución analítica vs BVProblem para la fuerza axial",
                        titlefont=font(10,"Computer Modern"),
                        shape = :star5, color = :blue )
                        

#Imprimo las dos graficos, en un subplot (layout (2,1))
plot(fig_despla, fig_axial, layout=(2,1))

#%%Fin
