# Elementos_Finitos
Códigos de análisis estructural por elementos finitos, traducidos de PYTHON y MATLAB,                                                                                              
(ver: https://github.com/diegoandresalvarez/elementosfinitos) al lenguaje de programación JULIA (https://julialang.org/). 

El archivo `Install_Pkg.jl` [clic](https://github.com/Sbeltranj/elementos-finitos/blob/master/Install_Pkg.jl), instala los paquetes de JULIA utilizados en este repositorio. 

Aquí encontrará algunos códigos, en los cuales se hace uso la librería de PYTHON (Matplotlib); ya que JULIA permite utilizar diversas librerías de diferentes lenguajes de programación, en su entorno de salida gráfica. (visitar : https://docs.juliaplots.org/latest/backends/).

Para ello se hace necesario, tener instalada la librería Matplotlib de Python: `pip install matplotlib`: 
https://matplotlib.org/stable/.

Además, de instalar los paquetes desde la consola de JULIA, para hacer los llamados a estos backend: `import Pkg` `Pkg.add("PyCall")`, `Pkg.add("PyPlot")`.
https://github.com/JuliaPy/PyPlot.jl#readme.

- Gancho analizado con EF isoparamétricos 8 nodos. [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/2D/Q8_Ejemplo)

[![git.png](https://i.postimg.cc/yxbTQTBs/git.png)](https://postimg.cc/p5K8yzp6)

- Momentos losa simplemente apoyada  [clic](https://github.com/Sbeltranj/elementos-finitos/tree/master/Losas/Kirchhoff_Love/Kirchhoff_Love).

[![plot-9.png](https://i.postimg.cc/m2vWxVsP/plot-9.png)](https://postimg.cc/vctNfrDy)
