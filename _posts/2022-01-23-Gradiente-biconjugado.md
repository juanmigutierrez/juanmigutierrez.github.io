---
layout: post
read_time: true
show_date: true
title:  Gradiente Biconjugado
date:   2022-01-23
description: Gradiente Biconjugado.
img: posts/20220323/gradiente.png
tags: [gradiente,biconjugado]
author: Miguel Gutierrez y David Felipe Martinez
Mathjax: yes
---
# Gradiente Conjugado

Intento 1

El método del gradiente biconjugado se basa en el método del gradiente conjugado donde se quiere solucionar el sistema <span class="math display">\[Ax = b\]</span>. Este método calcula un pseudo-gradiente <span class="math inline">\(\hat r_k\)</span> y una pseudo-dirección de descenso <span class="math inline">\(\hat p_k\)</span>. Estos pseudo gradientes serán ortogonales a los gradientes <span class="math inline">\(r_k\)</span> y <span class="math inline">\(p_k\)</span>. Este método a diferencia del gradiente conjugado no garantiza convergencia en <span class="math inline">\(m\)</span> pasos.  
El pre-condicionamiento mejora en general la convergencia de los métodos, esto pues transforma la matriz de coeficientes en una con mejor espectro.

**<span style="color: NavyBlue">Método de gradiente biconjugado precondicionado</span>**  
ENTRADA : el numero de ecuaciones y valores desconocidos <span class="math inline">\(n;\)</span> las entradas <span class="math inline">\(a_{ij}, 1\leq i,j\leq n\)</span> de la matriz <span class="math inline">\(A\)</span>; las entradas <span class="math inline">\(b_j,1\leq j \leq n\)</span> del vector <span class="math inline">\(\mathbf{b}\)</span>; las entradas <span class="math inline">\(\gamma_{ij}, 1\leq i,j\leq n\)</span> de la matriz precondicionada <span class="math inline">\(c^{-1}\)</span>, las entradas <span class="math inline">\(x_i,1\leq i\leq n\)</span> de la aproximación inicial <span class="math inline">\(\mathbf{x=x^{(0)}}\)</span>, el numero máximo de iteraciones <span class="math inline">\(N\)</span>, la tolerancia <span class="math inline">\(TOL\)</span>.  
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.  
_Paso 1_ Determine <span class="math inline">\(\mathbf{r=b-Ax;}\)</span> (_Calcule_ <span class="math inline">\(\mathbf{r^{(0)}}\)</span>)

<span class="math display">\[\begin{split} \mathbf{\hat r} &\mathbf{= r} \text{; (\textit{Calcule} $\mathbf{\hat r^{(0)}}$)} \\ \mathbf{q} &\mathbf{=C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(0)}}$)}\\ \mathbf{\hat q} &\mathbf{=(C^T)^{-1}\hat r}\text{; (\textit{Nota} $\mathbf{\hat q= \hat q^{(0)}}$)} \\ \end{split}\]</span>

_Paso 2_ Determine <span class="math inline">\(k=1\)</span>

_Paso 3_ Mientras (<span class="math inline">\(k \leq N\)</span>) haga los pasos 4-7

_Paso 4_ Si <span class="math inline">\(||r_k|| \leq TOL\)</span>, entonces, retornar.

_Paso 5_ Determine <span class="math display">\[\begin{split} \alpha &= \mathbf{\frac{\hat r \cdot q}{\hat p\cdot A \cdot p}} \text{; (\textit{Nota} $\mathbf{\hat \alpha= \frac{ \langle \hat r^{(k)}, q^{(k)}\rangle} {\langle \hat p^{(k)},A p^{(k)} \rangle}}$)}\\ \mathbf{r} &\mathbf{=r-\alpha q} \text{; (\textit{Nota} $\mathbf{r=r^{(k)}}$)}\\ \hat \alpha &= \mathbf{\frac{\hat r \cdot \hat q}{\hat p\cdot A^T \cdot \hat p}} \text{; (\textit{Nota} $\mathbf{\hat \alpha= \frac{ \langle \hat r^{(k)}, \hat q^{(k)}\rangle} {\langle \hat p^{(k)},A^T \hat p^{(k)} \rangle}}$)}\\ \mathbf{\hat r} & \mathbf{=\hat r-\hat \alpha \hat q} \\ \mathbf{x_{k+1}} &\mathbf{= x_k + \alpha p} \\ \end{split}\]</span>

_Paso 6_ Determine <span class="math display">\[\begin{split} \mathbf{q} & \mathbf{= C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(k)}}$)}\\ \mathbf{\hat q} &\mathbf{= (C^T)^{-1}\hat r} \text{; (\textit{Nota} $\mathbf{\hat q=\hat q^{(k)}}$)} \\ \mathbf{\beta} &\mathbf{= \frac{\hat r \cdot q}{\hat r \cdot q_{k-1}}} \text{; (\textit{Nota} $\mathbf{\beta= \frac{ \langle \hat r^{(k)}, q^{(k+1)}\rangle} {\langle \hat r^{(k)}, q^{(k)} \rangle}}$)} \end{split}\]</span>

_Paso 7_ Determine <span class="math display">\[\begin{split} \mathbf{p} &\mathbf{= q_{k+1}+\beta p} \text{; (\textit{Nota} $\mathbf{p=p^{(k)}}$)}\\ \mathbf{\hat p} &\mathbf{= \hat q_{k+1}+\beta \hat p} \text{; (\textit{Nota} $\mathbf{\hat p= \hat p^{(k)}}$)} \\ k &= k+1 \end{split}\]</span>

_Paso 8_ Si <span class="math inline">\((k>n)\)</span> entonces SALIDA. PARE.  
La secuencia de vectores satisfacen la condición de biortogonalidad

<span class="math display">\[\mathbf{\hat r_i \cdot r_j = r_i \cdot \hat r_j}= 0, \quad j<i\]</span>

y la condición de biconjugancia

<span class="math display">\[\mathbf{\hat p_i \cdot A \cdot p_j = p_i \cdot A^T \cdot \hat p_j}= 0, \quad j<i\]</span>

también se cumple ortogonalidad mutua

<span class="math display">\[\mathbf{\hat r_i \cdot p_j = r_i \cdot \hat p_j}= 0, \quad j<i\]</span>

Las demostraciones de estas propiedades resultan del proceso de inducción. Note que cuando <span class="math inline">\(A\)</span> es simétrica y se escoge <span class="math inline">\(\hat r = r\)</span>, luego <span class="math inline">\(\hat r_k = r_k\)</span> y <span class="math inline">\(\hat p_k = p_k\)</span> para todo <span class="math inline">\(k\)</span> es el caso del algoritmo del gradiente precondicionado. La matriz condicionadora mejora la convergencia del algoritmo <span class="citation" data-cites="Bi1"></span>. Note que el gradiente conjugado es el doble de costoso que el método del gradiente.

# Residuo Mínimo Generalizado (GMRES)

El método del Residuo Mínimo Generalizado (GMRES) es un método iterativo utilizado para encontrar la soluciona a un sistema lineal no simétrico. El método reside en reconstruir una base ortonormal del espacio de Krylov y por lo tanto es vulnerable a sesgos de errores computacionales. Este método reduce las dimensiones y por lo tanto resulta eficiente comparado con otros métodos.

El espacio de Krylov de dimensión <span class="math inline">\(k\)</span> se define como el subespacio lineal mas pequeño que contiene el conjunto por las imágenes de <span class="math inline">\(b\)</span> bajo las primeras <span class="math inline">\(k\)</span> potencias de A.

<span class="math display">\[\begin{split} \kappa_n &= span\{b,Ab,...,A^{n-1}b\}\\ &=\{q_1,q2,..,q_n\} \subset \mathbf{R}^{m} \end{split}\]</span>

Donde cada <span class="math inline">\(q_i\)</span> son vectores ortonormales. Ahora definimos la matriz de Krylov:

<span class="math display">\[K_n = \begin{bmatrix} \vdots & \vdots &\vdots &\vdots\\ b & Ab & A^2b & A^{n-1}b \\ \vdots & \vdots &\vdots &\vdots \end{bmatrix} = V_n R_n\]</span> Donde <span class="math inline">\(V_n\)</span> y <span class="math inline">\(R_n\)</span> conforman las matrices de una descomposición <span class="math inline">\(Q R\)</span> respectivamente.  
Cualquier vector <span class="math inline">\(x\in K_n\)</span> puede ser escrito como <span class="math inline">\(x = x_0 + V_n y_n\)</span> donde <span class="math inline">\(y\)</span> es un <span class="math inline">\(n\)</span>-vector y <span class="math inline">\(V_n\)</span> es una base ortonormal del subespacio de Krylov <span class="math inline">\(K_n\)</span>. Definimos la función de perdida que queremos minimizar.

<span class="math display">\[\begin{split} J(y) = ||b-Ax||_2 &=||b-A(x_0+V_n y)||_2\\ &= ||r^0-AV_n y||_2 \\ \end{split}\]</span>

Con la relación de arnoldi <span class="math display">\[AV_n = V_{n+1}\bar H_k\]</span>

obtenemos

<span class="math display">\[\begin{split} \left\|r_{0}-A V_{n} y\right\|_{2}&=\left\|r_{0}-V_{n+1} \bar{H}_{n} y\right\|_{2}\\ &=\left\|V_{n+1}\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\\ &=\left\|\beta e_{1}-\bar{H}_{n} y\right\|_{2} \end{split}\]</span>

Donde <span class="math inline">\(r_0 = \beta v_1\)</span>, <span class="math inline">\(e_1\)</span> es la primera columna de la matriz identidad <span class="math inline">\((k+1)\times(k+1)\)</span> y <span class="math inline">\(\bar H\)</span> es una matriz superior de Hessenberg

<span class="math display">\[\bar H_{n}=\left[\begin{array}{cccc} h_{11} & h_{12} & \cdots & h_{1 n} \\ h_{21} & h_{22} & \cdots & h_{2 n} \\ 0 & h_{32} & \ddots & h_{3 n} \\ \vdots & \vdots & & \vdots \\ 0 & 0 & \cdots & h_{n n} \end{array}\right]\]</span>

Luego el problema se convierte en minimizar en cada iteración

<span class="math display">\[\begin{split} x_{n}&=x_{0}+V_{n} y_{n} \\ y_{n}&=\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2} \end{split}\]</span>

Para encontrar <span class="math inline">\(V_n\)</span> el algoritmo de Arnoldi nos será de utilidad.  
**<span style="color: NavyBlue">Método de residuo Mínimo Generalizado (GMRES)</span>**  
ENTRADA : <span class="math inline">\(x_0,b,A\)</span>  
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.  
_Paso 1_ Determine <span class="math inline">\(\mathbf{r_0=b-Ax_0}\)</span> y <span class="math inline">\(v_1 = \mathbf{r_0/||r_0||_2}\)</span>

_Paso 2_ Defina una matriz superior de Hesenberg <span class="math inline">\(H_n=0\)</span> de tamaño <span class="math inline">\((m+1)\times m\)</span>.

_Paso 3_ Para <span class="math inline">\(j=1,2,...,m\)</span> haga los pasos 4, 5, 6 y 7

_Paso 4_ Determine <span class="math inline">\(w_j = Av_j\)</span>

_Paso 5_ Para <span class="math inline">\(i=1,2,...,m\)</span> determine <span class="math display">\[\begin{split} h_{ij} &= \langle w_j,vi \rangle \\ w_j &= w_j - h_{ij}v_i \end{split}\]</span>

_Paso 6_ Determine <span class="math inline">\(h_{j+1j} = ||w_j||_2\)</span>.

_Paso 7_ Si <span class="math inline">\(h_{j+1,j}=0\)</span> entonces <span class="math inline">\(m=j\)</span> y romper iteración (break).

_Paso 8_ Determine <span class="math inline">\(v_{j+1} = w_j/h_{j+1,j}\)</span>.

_Paso 9_ Computar minimizador de <span class="math inline">\(\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\)</span> y <span class="math inline">\(x_n=x_0+V_ny_n\)</span>.  
Note que en los pasos <span class="math inline">\(1-8\)</span> se esta encontrando una base ortonormal del espacio de Krylov. De aquí en adelante, la pregunta ¿Como encontramos <span class="math inline">\(y_m\)</span> ? Es donde surgen diferentes métodos. Uno de los mas simple y poderosos fue introducido por Saad Y Schultz.

## GMRES Estandar

Para solucionar el problema de mínimos cuadrados

<span class="math display">\[y_{n} =\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\]</span>

Se sugirió utilizar unas rotaciones dadas para convertir la anterior ecuación en un sistema superior triangular para descomponer <span class="math inline">\(H_n\)</span> en <span class="math inline">\(Q_n \hat R_k\)</span> Donde <span class="math inline">\(Q_n\)</span> es una matriz ortonormal obtenida del producto de <span class="math inline">\(k\)</span> rotaciones y <span class="math inline">\(\hat R_k =\begin{bmatrix}R_k \\ 0...0 \end{bmatrix}\)</span>. La matriz de rotación es una matriz identidad de orden <span class="math inline">\(n+1\)</span> donde solo cuatro componentes son remplazadas por escalares <span class="math inline">\(c_i,s_i\)</span> de la siguiente manera

<span class="math display">\[J_i = \begin{bmatrix} 1 & & & 0 \\ &c_i & s_i & \\ & -s_i &c_i & \\ 0 & & & 1 \\ \end{bmatrix} \quad \text{donde } c_i^2+s_i^2 = 1,i= 1,...,n\]</span>

La matriz es construida de tal manera que

<span class="math display">\[\begin{bmatrix} c_i & s_i \\ -s_i &c_i \\ \end{bmatrix} \begin{bmatrix} h_{i,i} \\ h_{i+1,i} \\ \end{bmatrix} = \begin{bmatrix} * \\ 0 \\ \end{bmatrix}\]</span>

El producto punto de las primeras <span class="math inline">\(n\)</span> rotaciones, es decir, <span class="math inline">\(Q_n = J_nJ_{k-1}\hdots J_1\)</span> al aplicarse sobre la matriz por izquierda de <span class="math inline">\(H_n\)</span> y <span class="math inline">\(\beta e_1\)</span> donde se obtiene la descomposición

<span class="math display">\[R_ny_n = g_n\]</span>

y la norma residual se obtiene

<span class="math display">\[\begin{split} \left(\bar{H}_{k} \beta e_{1}\right) =&\left(\begin{array}{cccccc} * & * & * & * & * & \beta \\ * & * & * & * & * & 0 \\ & * & * & * & * & 0 \\ & & * & * & * & 0 \\ & & & * & * & 0 \\ & & & & * & 0 \end{array}\right) \\ &\frac{\text { Aplicando las }}{\text {rotaciones } \overrightarrow{Q_{k}}} \\ &\left(\begin{array}{cccccc} + & + & + & + & + & \times \\ 0 & + & + & + & + & \times \\ & 0 & + & + & + & \times \\ & & 0 & + & + & \times \\ & & & 0 & + & \times \\ & & & & 0 & \gamma_{k} \end{array}\right)=\left(\begin{array}{cc} R_{n} & g_{n} \\ 0 & \gamma_{n} \end{array}\right) \end{split}\]</span>

Sujeto a que <span class="math inline">\(g_n^T \in R^n\)</span>. El algoritmo finaliza cuando <span class="math inline">\(\gamma_k < \varepsilon\)</span>. <span class="citation" data-cites="Bi3"></span><span class="citation" data-cites="Bi4"></span><span class="citation" data-cites="Bi5"></span>. Finalmente, el algoritmo después de hallar los vector ortonormales con Arnoldi quedaría:  
**<span style="color: NavyBlue">Método Estándar (GMRES)</span>**  

ENTRADA : <span class="math inline">\(x_0,bA\)</span>  
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.  
_Paso 1_ Determine <span class="math inline">\(\mathbf{r_0=b-A x_0}\)</span> y <span class="math inline">\(v_1 = \mathbf{r_0/||r_0||_2}\)</span>

_Paso 2_ Para <span class="math inline">\(i=1,2,...,m\)</span> haga los pasos 3, 4, 5, y 6

_Paso 3_ <span class="math inline">\(v_{j+1} = \prod_j A v_j / ||A v_j||\)</span>

_Paso 4_ <span class="math inline">\(\begin{bmatrix}R_j \\ 0 \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1\bar H_j)\)</span>

_Paso 5_ <span class="math inline">\(\begin{bmatrix}g_j \\ \gamma_j \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1(\beta e_1))\)</span>

_Paso 6_ Si <span class="math inline">\(\gamma_j < T O L\)</span>, break y haga el paso 7 y 8

_Paso 7_ <span class="math inline">\(y_k = R_k^{-1}g_k\)</span> y <span class="math inline">\(x_k = x_0 + v_k y_k\)</span>

_Paso 8_ Si <span class="math inline">\(x_k\)</span> no es solución, entonces determine <span class="math inline">\(x_0=x_k\)</span> y retome el paso 1.  
Observe que el paso 1, y 3 esta mejor explicado en el primer pseudo-código, por lo que este ultimo se centra en la ultima parte del algoritmo que es encontrar <span class="math inline">\(y_m\)</span> y finalmente obtener <span class="math inline">\(x\)</span>.  
Además, podemos obtener las siguientes afirmaciones para los GMRES estándar [2, teorema 4]:

*   El rango de <span class="math inline">\(A V_k\)</span> es igual al rango de <span class="math inline">\(R_k\)</span>, en particular, si <span class="math inline">\(r_{k,k}=0\)</span>, entonces A es singular.

*   El vector <span class="math inline">\(y_k\)</span> que minimiza <span class="math inline">\(\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\)</span> y <span class="math inline">\(x_n=x_0+V_ny_n\)</span>, esta dado por <span class="math inline">\(y_k=R_k^{-1}g_k\)</span>

*   la normal del residual en el paso <span class="math inline">\(k\)</span> satisface <span class="math inline">\(||b - A x_k|| = |\gamma_k|\)</span>

# Complejidades

La complejidad computacional de GMRES <span class="math inline">\(O(n^2+n k+k^2)\)</span> y de memoria es <span class="math inline">\(O(n^2+n k+k^2)\)</span> donde el <span class="math inline">\(k\)</span> el numero de pasos de ortogonalización de. Mientras lo que se conoce del algoritmo Biconjugado es <span class="math inline">\(O(n)\)</span> en cada iteración.