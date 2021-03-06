using Interpolations, PyPlot, PyCall, MappedArrays, CubicSplines

# https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/1D/EF_barra_n_nodos/comparing_interpolation_algorithms.m
#https://juliamath.github.io/Interpolations.jl/latest/api/

ENV["MPLBACKEND"]="qt5agg"
pygui(true)

## Se crea un espacio para hacer clic y definir los puntos de interpolación 
figure()
title("Haciendo clic con el ratón, defina los 10 puntos a interpolar")
axis([-6, 6, -6, 6])

coord = ginput(10)
xi   = collect(1:1:10)
xnod = mappedarray(coord->first(coord),coord)
ynod = mappedarray(coord->last(coord),coord)


itp = interpolate(ynod, BSpline(Linear()))
#itp1 = interpolate(ynod, BSpline(NoInterp()))
#itp = CubicSplineInterpolation(ynod, xi)

spline = CubicSpline(xnod, ynod)
xs = range(xnod[1], stop=xnod[end], length=1000)
ys = spline[xs]

plt.grid()
plot(xnod, ynod, "o")
plot(xnod, itp)
plot(xs, ys)
#plot(xnod, itp1)

##SIN TERMINAR-FALTA FALTA

#= using Dierckx
fi = CubicSplineInterpolation(xnod, ynod, extrapolation_bc=Flat())  # use Interpolations
fd = Spline1D(xs, ys, bc="nearest")  # use Dierckx
spline = Spline1D(xnod,ynod,k=3,bc = "nearest" =#)