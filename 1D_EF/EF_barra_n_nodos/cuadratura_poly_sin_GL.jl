# Programa elaborado en JULIA 1.7.1

# Diego Andrés Alvarez Marín
# daalvarez@unal.edu.co
# https://github.com/diegoandresalvarez/elementosfinitos/tree/master/codigo/1D/EF_barra_n_nodos

# Traducido por:
# Santiago Beltrán Jaramillo
# sbeltran@unal.edu.co

using Plots
using Polynomials

# realiza integraciones de funciones utilizando cuadraturas de Gauss-Legendre
# comparando los errores obtenidos con la solución exacta.
include("gausslegendre_quad.jl")

plotlyjs()   #Pkg.add("PlotlyJS")

# Integración de:
#
f(x)  = 0.2 + 25*x - 200*x.^2 +675*x.^3 - 900*x.^4 + 400*x.^5
#
# entre 0 y 0.8 usando cuadraturas de Gauss-Legendre

a   = 0;                  # límites de integración
b   = 0.8;        
sol = 3076/1875;          # solución exacta
err = zeros(10,1);        # separo la memoria
for m = 1:10              # varío el número de puntos de la cuadratura
    w_gl, x_gl = gausslegendre_quad(m);  # calculo w y c de la cuadratura
     err[m] = abs(((b-a)/2)*sum(w_gl.*f.((b+a)/2 .+ (b-a)*x_gl/2)) - sol);
end

n_1 = plot()
n_1 = plot(err, title = "Cuadratura de Gauss Legendre",
           yaxis = "Error", xaxis = "Número de puntos en la cuadratura",
           label = nothing )

# Integracion de:
#
f(x)   = sin(x)
#
# entre 0 y pi/2 usando cuadraturas de Gauss-Legendre

a   = 0;                  # límites de integracion
b   = pi/2;
sol = 1;                  # solución exacta
err = zeros(10,1);        # separo la memoria
for m = 1:10
   w_gl, x_gl = gausslegendre_quad(m);
   err[m] = abs(((b-a)/2)*sum(w_gl.*f.((b+a)/2 .+ (b-a)*x_gl/2)) - sol);
end

n_2 = plot()

n_2 = plot(err, yaxis=:log, title = "Cuadratura de Gauss Legendre",
           xaxis = "Número de puntos en la cuadratura", label = nothing )
ylabel!("Error")
xlabel!("Número de puntos en la cuadratura")

display(n_1)
display(n_2)




