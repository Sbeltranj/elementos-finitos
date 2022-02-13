# Elementos_Finitos
Códigos de análisis estructural por elementos finitos, traducidos de PYTHON y MATLAB,                                                                                              
(ver: https://github.com/diegoandresalvarez/elementosfinitos) al lenguaje de programación JULIA (https://julialang.org/). 

El archivo `Install_Pkg.jl` [clic](https://github.com/Sbeltranj/elementos-finitos/blob/master/Install_Pkg.jl), instala los paquetes de JULIA utilizados en este repositorio. 

## Contenido 

- Elementos finitos de barra (1D).  [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/1D_EF)
- Análisis de estructuras en tensión plana, con EFs T3, Q8, QM6(empleado en SAP2000), Q6. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/2D)
- Método matricial (Pórticos, cerchas, vigas). [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/An%C3%A1lisis%20matricial)
- Elementos finitos tridimensionales H20. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/3D)
- Losas macizas, por las teorías de Kirchhoff-Love y Mindlin–Reissner. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/Losas)
- Elementos finitos de viga, por las teorías de Euler-Bernoulli y Timoshenko-Ehrenfest. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/Vigas)

## PyPlot Backend 

Aquí encontrará algunos códigos, en los cuales se hace uso la librería de PYTHON (Matplotlib); ya que JULIA permite utilizar diversas librerías de diferentes lenguajes de programación, en su entorno de salida gráfica. (visitar : https://docs.juliaplots.org/latest/backends/).

Para ello se hace necesario, tener instalada la librería Matplotlib de Python: `pip install matplotlib`: 
https://matplotlib.org/stable/.

Además, de instalar los paquetes desde la consola de JULIA, para hacer los llamados a estos backend: `import Pkg` `Pkg.add("PyCall")`, `Pkg.add("PyPlot")`.
https://github.com/JuliaPy/PyPlot.jl#readme.

- Gancho analizado con EF isoparamétricos 8 nodos. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/2D/Q8_Ejemplo)

[![git.png](https://i.postimg.cc/yxbTQTBs/git.png)](https://postimg.cc/p5K8yzp6)

- Momentos flectores y torsores, losa maciza simplemente apoyada [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/Losas/Kirchhoff_Love/Kirchhoff_Love).

[![plot-9.png](https://i.postimg.cc/m2vWxVsP/plot-9.png)](https://postimg.cc/vctNfrDy)

- Modos de energía nula, EF QL9.[clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/Losas/MINDLIN/QL9_integracion_reducida).

[![MEN.png](https://i.postimg.cc/sgBsG3FL/MEN.png)](https://postimg.cc/qNd99Vvc)
