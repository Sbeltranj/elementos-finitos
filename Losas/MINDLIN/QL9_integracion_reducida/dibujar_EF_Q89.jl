function dibujar_EF_Q89_RM(xe, ye, N, ae, te, esc_w, esc_t)
    
    # Programa original (MATLAB) elaborado por:
    # Diego Andrés Alvarez Marín
    # daalvarez@unal.edu.co

    #Traduciendo a JULIA 1.7.1, por:
    # Santiago Beltrán Jaramillo

    #Importante: este programa hace uso de PyPlot de Python, importante tener dicha librería instalada.

    ## Dibuja un EF de losa de Mindlin de 8 o 9 nodos QL9 deformado
    # xe, ye  nodos del borde de la losa
    # N       funciones de forma de la losa
    # ae      desplazamientos nodales en los nodos del borde de la losa
    # te      espesor de la losa
    # esc_w   escalamiento del desplazamiento vertical
    # esc_t   escalamiento de los giros

    ## Se definen algunas constantes
    w_ = 1;  tx_ = 2;  ty_ = 3; # lectura del código
    
    nno = size(xe,1);           # número de nodos del elemento finito
    
    ## Se calcula la geometría isoparamétrica
    xi  = [       collect(-1:0.1:1)'  +ones(1,length(-1:0.1:1))                collect(1:-0.1:-1)'   -ones(1,length(1:-0.1:-1)) ];
    eta = [-ones(1,length(-1:0.1:1))               collect(-1:0.1:1)'    +ones(1,length(1:-0.1:-1))                collect(1:-0.1:-1)'   ];
    
    n = length(xi);  # number of points of vector xi
    
    x = zeros(n,1);
    y = zeros(n,1);
    for j = 1:n
       x[j] = sum(N(xi[j], eta[j]) .* xe);
       y[j] = sum(N(xi[j], eta[j]) .* ye);
    end
    
    ## Se gráfica el elemento original


    ax.plot3D(vec(x), vec(y), repeat([-te/2],size(x)[1]), "b-")
    ax.plot3D(vec(x), vec(y), zeros(size(x)[1]), "b-")
    ax.plot3D(vec(x), vec(y), repeat([te/2],size(x)[1]), "b-") 
    
    ## Se grafican las lineas verticales originales
    for j = 1:nno
        ax.plot3D([ xe[j], xe[j] ], [ ye[j], ye[j] ], [ -te/2, te/2 ], "b*-");
    end 
    
    ## Se calcula el elemento deformado
    mov = reshape(ae,3,nno)';  # se reorganiza el vector de movimientos nodales

    mov[:,w_]        = esc_w*mov[:,w_];        # se escala la deform vertical
    mov[:,[tx_ ty_]] = esc_t*mov[:,[tx_ ty_]]; # se escalan los ángulos
    
    w  = zeros(n,1);
    tx = zeros(n,1);
    ty = zeros(n,1);

    for j = 1:n
       w[j]  = sum(N(xi[j], eta[j]) .* mov[:,w_]);
       tx[j] = sum(N(xi[j], eta[j]) .* mov[:,tx_]);
       ty[j] = sum(N(xi[j], eta[j]) .* mov[:,ty_]);   
    end
    z =  te/2; u_sup = -z*tx; v_sup = -z*ty; w_sup = w.+z;
    z =     0; u_mid = -z*tx; v_mid = -z*ty; w_mid = w.+z;
    z = -te/2; u_inf = -z*tx; v_inf = -z*ty; w_inf = w.+z;
    
    ## Se gráfica el elemento deformado
    ax.plot3D(vec(x+u_sup), vec(y+v_sup), vec(w_sup), "r"); # cara superior
    ax.plot3D(vec(x+u_mid), vec(y+v_mid), vec(w_mid), "r"); # plano medio
    ax.plot3D(vec(x+u_inf), vec(y+v_inf), vec(w_inf), "r"); # cara inferior 
    
    ## Se grafican las lineas verticales deformadas
    z =  te/2; u_sup = -z*mov[:,tx_]; v_sup = -z*mov[:,ty_]; w_sup = mov[:,w_].+z;
    z = -te/2; u_inf = -z*mov[:,tx_]; v_inf = -z*mov[:,ty_]; w_inf = mov[:,w_].+z;
    
    for j = 1:nno
        ax.plot3D([ xe[j]+u_inf[j], xe[j]+u_sup[j] ], [ ye[j]+v_inf[j], ye[j]+v_sup[j] ], [ w_inf[j], w_sup[j] ], "r*-");
    end 
    
    ## Meshgrid, tomado de: https://stackoverflow.com/questions/44581049/utilizing-ndgrid-meshgrid-functionality-in-julia
    function meshgrid(x, y)
        X = [i for i in x, j in 1:length(y)]
        Y = [j for i in 1:length(x), j in y]
        return X, Y
    end

    delta = 0.1
    xxi  = collect(-1:delta:1)
    eeta = collect(-1:delta:1)
    n    = length(xxi);
    xi, eta = meshgrid(xxi, eeta)
    xi = xi' ; eta = eta'
    

    w = zeros(21,21);
    x = zeros(21,21);
    y = zeros(21,21);
    for i = 1:21
       for j = 1:21
          w[i,j] = sum(N(xi[i,j], eta[i,j]) .* mov[:,w_]);
          x[i,j] = sum(N(xi[i,j], eta[i,j]) .* xe);
          y[i,j] = sum(N(xi[i,j], eta[i,j]) .* ye);
       end
    end
    
    #@pyimport matplotlib.cm as cmapi
    cm = pyimport("matplotlib.cm")

#=     norma = plt.Normalize(minimum(w), maximum(w))
    cim = cm.jet(norma(w)) 
    rcount, ccount, _ = size(cim) =#

    #my_col = cm.jet(w/minimum(w))
    ax.plot_surface(x, y, w, rstride=1, cstride=1, 
        linewidth=2, cmap = "jet")

    return
end
