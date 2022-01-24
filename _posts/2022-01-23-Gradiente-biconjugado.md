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
mathjax: yes
---

::: {#title-block-header}
# Métodos del Gradiente Biconjugado y Mínimos Cuadrados Generalizados GMRES

# Gradiente Conjugado

<p>El método del gradiente biconjugado se basa en el método del gradiente conjugado donde se quiere solucionar el sistema <span class="math display">\[Ax = b\]</span>. Este método calcula un pseudo-gradiente <span class="math inline">\(\hat r_k\)</span> y una pseudo-dirección de descenso <span class="math inline">\(\hat p_k\)</span>. Estos pseudo gradientes serán ortogonales a los gradientes <span class="math inline">\(r_k\)</span> y <span class="math inline">\(p_k\)</span>. Este método a diferencia del gradiente conjugado no garantiza convergencia en <span class="math inline">\(m\)</span> pasos.<br />
El pre-condicionamiento mejora en general la convergencia de los métodos, esto pues transforma la matriz de coeficientes en una con mejor espectro.</p>
<p><strong><span style="color: NavyBlue">Método de gradiente biconjugado precondicionado</span></strong><br />
ENTRADA : el numero de ecuaciones y valores desconocidos <span class="math inline">\(n;\)</span> las entradas <span class="math inline">\(a_{ij}, 1\leq i,j\leq n\)</span> de la matriz <span class="math inline">\(A\)</span>; las entradas <span class="math inline">\(b_j,1\leq j \leq n\)</span> del vector <span class="math inline">\(\mathbf{b}\)</span>; las entradas <span class="math inline">\(\gamma_{ij}, 1\leq i,j\leq n\)</span> de la matriz precondicionada <span class="math inline">\(c^{-1}\)</span>, las entradas <span class="math inline">\(x_i,1\leq i\leq n\)</span> de la aproximación inicial <span class="math inline">\(\mathbf{x=x^{(0)}}\)</span>, el numero máximo de iteraciones <span class="math inline">\(N\)</span>, la tolerancia <span class="math inline">\(TOL\)</span>.<br />
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.<br />
<em>Paso 1</em> Determine <span class="math inline">\(\mathbf{r=b-Ax;}\)</span> (<em>Calcule</em> <span class="math inline">\(\mathbf{r^{(0)}}\)</span>)</p>
<p><span class="math display">\[\begin{split}
   \mathbf{\hat r} &amp;\mathbf{= r} \text{; (\textit{Calcule} $\mathbf{\hat r^{(0)}}$)} \\
   \mathbf{q} &amp;\mathbf{=C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(0)}}$)}\\
   \mathbf{\hat q} &amp;\mathbf{=(C^T)^{-1}\hat r}\text{; (\textit{Nota} $\mathbf{\hat q= \hat q^{(0)}}$)} \\
\end{split}\]</span></p>
<p><em>Paso 2</em> Determine <span class="math inline">\(k=1\)</span></p>
<p><em>Paso 3</em> Mientras (<span class="math inline">\(k \leq N\)</span>) haga los pasos 4-7</p>
<p><em>Paso 4</em> Si <span class="math inline">\(||r_k|| \leq TOL\)</span>, entonces, retornar.</p>
<p><em>Paso 5</em> Determine <span class="math display">\[\begin{split}
   \alpha &amp;= \mathbf{\frac{\hat r \cdot q}{\hat p\cdot A \cdot p}}  \text{; (\textit{Nota} $\mathbf{\hat \alpha=  \frac{ \langle \hat r^{(k)},  q^{(k)}\rangle} {\langle \hat p^{(k)},A p^{(k)} \rangle}}$)}\\
   \mathbf{r} &amp;\mathbf{=r-\alpha q} \text{; (\textit{Nota} $\mathbf{r=r^{(k)}}$)}\\
    \hat \alpha &amp;= \mathbf{\frac{\hat r \cdot \hat q}{\hat p\cdot A^T \cdot \hat p}} \text{; (\textit{Nota} $\mathbf{\hat \alpha=  \frac{ \langle \hat r^{(k)}, \hat q^{(k)}\rangle} {\langle \hat p^{(k)},A^T \hat p^{(k)} \rangle}}$)}\\
   \mathbf{\hat r} &amp; \mathbf{=\hat r-\hat \alpha \hat q} \\
   \mathbf{x_{k+1}} &amp;\mathbf{= x_k + \alpha p} \\
\end{split}\]</span></p>
<p><br />
<em>Paso 6</em> Determine <span class="math display">\[\begin{split}
   \mathbf{q} &amp; \mathbf{= C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(k)}}$)}\\
   \mathbf{\hat q} &amp;\mathbf{= (C^T)^{-1}\hat r} \text{; (\textit{Nota} $\mathbf{\hat q=\hat q^{(k)}}$)} \\
   \mathbf{\beta} &amp;\mathbf{= \frac{\hat r \cdot q}{\hat r \cdot q_{k-1}}} \text{; (\textit{Nota} $\mathbf{\beta=  \frac{ \langle \hat r^{(k)}, q^{(k+1)}\rangle} {\langle \hat r^{(k)},   q^{(k)} \rangle}}$)}
\end{split}\]</span></p>
<p><em>Paso 7</em> Determine <span class="math display">\[\begin{split}
   \mathbf{p} &amp;\mathbf{= q_{k+1}+\beta p} \text{; (\textit{Nota} $\mathbf{p=p^{(k)}}$)}\\
   \mathbf{\hat p} &amp;\mathbf{= \hat q_{k+1}+\beta \hat p} \text{; (\textit{Nota} $\mathbf{\hat p= \hat p^{(k)}}$)} \\
   k &amp;= k+1
\end{split}\]</span></p>
<p><em>Paso 8</em> Si <span class="math inline">\((k&gt;n)\)</span> entonces SALIDA. PARE.<br />
La secuencia de vectores satisfacen la condición de biortogonalidad</p>
<p><span class="math display">\[\mathbf{\hat r_i \cdot r_j = r_i \cdot \hat r_j}= 0, \quad j&lt;i\]</span></p>
<p>y la condición de biconjugancia</p>
<p><span class="math display">\[\mathbf{\hat p_i \cdot A \cdot p_j = p_i \cdot A^T \cdot \hat p_j}= 0, \quad j&lt;i\]</span></p>
<p>también se cumple ortogonalidad mutua</p>
<p><span class="math display">\[\mathbf{\hat r_i \cdot p_j = r_i \cdot \hat p_j}= 0, \quad j&lt;i\]</span></p>
<p>Las demostraciones de estas propiedades resultan del proceso de inducción. Note que cuando <span class="math inline">\(A\)</span> es simétrica y se escoge <span class="math inline">\(\hat r = r\)</span>, luego <span class="math inline">\(\hat r_k = r_k\)</span> y <span class="math inline">\(\hat p_k = p_k\)</span> para todo <span class="math inline">\(k\)</span> es el caso del algoritmo del gradiente precondicionado. La matriz condicionadora mejora la convergencia del algoritmo <span class="citation" data-cites="Bi1"></span>. Note que el gradiente conjugado es el doble de costoso que el método del gradiente.</p>
<h1 id="residuo-mínimo-generalizado-gmres">Residuo Mínimo Generalizado (GMRES)</h1>
<p>El método del Residuo Mínimo Generalizado (GMRES) es un método iterativo utilizado para encontrar la soluciona a un sistema lineal no simétrico. El método reside en reconstruir una base ortonormal del espacio de Krylov y por lo tanto es vulnerable a sesgos de errores computacionales. Este método reduce las dimensiones y por lo tanto resulta eficiente comparado con otros métodos.</p>
<p>El espacio de Krylov de dimensión <span class="math inline">\(k\)</span> se define como el subespacio lineal mas pequeño que contiene el conjunto por las imágenes de <span class="math inline">\(b\)</span> bajo las primeras <span class="math inline">\(k\)</span> potencias de A.</p>
<p><span class="math display">\[\begin{split}
    \kappa_n &amp;= span\{b,Ab,...,A^{n-1}b\}\\
    &amp;=\{q_1,q2,..,q_n\}  \subset \mathbf{R}^{m}
\end{split}\]</span></p>
<p>Donde cada <span class="math inline">\(q_i\)</span> son vectores ortonormales. Ahora definimos la matriz de Krylov:</p>
<p><span class="math display">\[K_n =
\begin{bmatrix}
\vdots &amp; \vdots &amp;\vdots &amp;\vdots\\
b &amp; Ab &amp; A^2b &amp; A^{n-1}b \\
\vdots &amp; \vdots &amp;\vdots &amp;\vdots
\end{bmatrix}
= V_n R_n\]</span> Donde <span class="math inline">\(V_n\)</span> y <span class="math inline">\(R_n\)</span> conforman las matrices de una descomposición <span class="math inline">\(Q R\)</span> respectivamente.<br />
Cualquier vector <span class="math inline">\(x\in K_n\)</span> puede ser escrito como <span class="math inline">\(x = x_0 + V_n y_n\)</span> donde <span class="math inline">\(y\)</span> es un <span class="math inline">\(n\)</span>-vector y <span class="math inline">\(V_n\)</span> es una base ortonormal del subespacio de Krylov <span class="math inline">\(K_n\)</span>. Definimos la función de perdida que queremos minimizar.</p>
<p><span class="math display">\[\begin{split}
    J(y) = ||b-Ax||_2 &amp;=||b-A(x_0+V_n y)||_2\\
    &amp;= ||r^0-AV_n y||_2 \\
\end{split}\]</span></p>
<p>Con la relación de arnoldi <span class="math display">\[AV_n = V_{n+1}\bar H_k\]</span></p>
<p>obtenemos</p>
<p><span class="math display">\[\begin{split}
\left\|r_{0}-A V_{n} y\right\|_{2}&amp;=\left\|r_{0}-V_{n+1} \bar{H}_{n} y\right\|_{2}\\
&amp;=\left\|V_{n+1}\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\\
&amp;=\left\|\beta e_{1}-\bar{H}_{n} y\right\|_{2}
\end{split}\]</span></p>
<p>Donde <span class="math inline">\(r_0 = \beta v_1\)</span>, <span class="math inline">\(e_1\)</span> es la primera columna de la matriz identidad <span class="math inline">\((k+1)\times(k+1)\)</span> y <span class="math inline">\(\bar H\)</span> es una matriz superior de Hessenberg</p>
<p><span class="math display">\[\bar H_{n}=\left[\begin{array}{cccc}
h_{11} &amp; h_{12} &amp; \cdots &amp; h_{1 n} \\
h_{21} &amp; h_{22} &amp; \cdots &amp; h_{2 n} \\
0 &amp; h_{32} &amp; \ddots &amp; h_{3 n} \\
\vdots &amp; \vdots &amp; &amp; \vdots \\
0 &amp; 0 &amp; \cdots &amp; h_{n n}
\end{array}\right]\]</span></p>
<p>Luego el problema se convierte en minimizar en cada iteración</p>
<p><span class="math display">\[\begin{split}
x_{n}&amp;=x_{0}+V_{n} y_{n} \\
y_{n}&amp;=\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}
\end{split}\]</span></p>
<p>Para encontrar <span class="math inline">\(V_n\)</span> el algoritmo de Arnoldi nos será de utilidad.<br />
<strong><span style="color: NavyBlue">Método de residuo Mínimo Generalizado (GMRES)</span></strong><br />
ENTRADA : <span class="math inline">\(x_0,b,A\)</span><br />
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.<br />
<em>Paso 1</em> Determine <span class="math inline">\(\mathbf{r_0=b-Ax_0}\)</span> y <span class="math inline">\(v_1 = \mathbf{r_0/||r_0||_2}\)</span></p>
<p><em>Paso 2</em> Defina una matriz superior de Hesenberg <span class="math inline">\(H_n=0\)</span> de tamaño <span class="math inline">\((m+1)\times m\)</span>.</p>
<p><em>Paso 3</em> Para <span class="math inline">\(j=1,2,...,m\)</span> haga los pasos 4, 5, 6 y 7</p>
<p><em>Paso 4</em> Determine <span class="math inline">\(w_j = Av_j\)</span></p>
<p><em>Paso 5</em> Para <span class="math inline">\(i=1,2,...,m\)</span> determine <span class="math display">\[\begin{split}
    h_{ij} &amp;= \langle w_j,vi \rangle \\
    w_j &amp;= w_j - h_{ij}v_i
\end{split}\]</span></p>
<p><em>Paso 6</em> Determine <span class="math inline">\(h_{j+1j} = ||w_j||_2\)</span>.</p>
<p><em>Paso 7</em> Si <span class="math inline">\(h_{j+1,j}=0\)</span> entonces <span class="math inline">\(m=j\)</span> y romper iteración (break).</p>
<p><em>Paso 8</em> Determine <span class="math inline">\(v_{j+1} = w_j/h_{j+1,j}\)</span>.</p>
<p><em>Paso 9</em> Computar minimizador de <span class="math inline">\(\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\)</span> y <span class="math inline">\(x_n=x_0+V_ny_n\)</span>.<br />
Note que en los pasos <span class="math inline">\(1-8\)</span> se esta encontrando una base ortonormal del espacio de Krylov. De aquí en adelante, la pregunta ¿Como encontramos <span class="math inline">\(y_m\)</span> ? Es donde surgen diferentes métodos. Uno de los mas simple y poderosos fue introducido por Saad Y Schultz.</p>
<h2 id="gmres-estandar">GMRES Estandar</h2>
<p>Para solucionar el problema de mínimos cuadrados</p>
<p><span class="math display">\[y_{n} =\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\]</span></p>
<p>Se sugirió utilizar unas rotaciones dadas para convertir la anterior ecuación en un sistema superior triangular para descomponer <span class="math inline">\(H_n\)</span> en <span class="math inline">\(Q_n \hat R_k\)</span> Donde <span class="math inline">\(Q_n\)</span> es una matriz ortonormal obtenida del producto de <span class="math inline">\(k\)</span> rotaciones y <span class="math inline">\(\hat R_k =\begin{bmatrix}R_k \\ 0...0
\end{bmatrix}\)</span>. La matriz de rotación es una matriz identidad de orden <span class="math inline">\(n+1\)</span> donde solo cuatro componentes son remplazadas por escalares <span class="math inline">\(c_i,s_i\)</span> de la siguiente manera</p>
<p><span class="math display">\[J_i =
\begin{bmatrix}
1 &amp; &amp; &amp; 0 \\
&amp;c_i &amp; s_i &amp; \\ 
&amp; -s_i &amp;c_i &amp; \\
0 &amp; &amp; &amp; 1 \\
\end{bmatrix}
\quad \text{donde } c_i^2+s_i^2 = 1,i= 1,...,n\]</span></p>
<p>La matriz es construida de tal manera que</p>
<p><span class="math display">\[\begin{bmatrix}
c_i &amp; s_i  \\ 
-s_i &amp;c_i \\
\end{bmatrix}
\begin{bmatrix}
h_{i,i}  \\ 
h_{i+1,i} \\
\end{bmatrix}
= 
\begin{bmatrix}
* \\ 
0 \\
\end{bmatrix}\]</span></p>
<p>El producto punto de las primeras <span class="math inline">\(n\)</span> rotaciones, es decir, <span class="math inline">\(Q_n = J_nJ_{k-1}\hdots J_1\)</span> al aplicarse sobre la matriz por izquierda de <span class="math inline">\(H_n\)</span> y <span class="math inline">\(\beta e_1\)</span> donde se obtiene la descomposición</p>
<p><span class="math display">\[R_ny_n = g_n\]</span></p>
<p>y la norma residual se obtiene</p>
<p><span class="math display">\[\begin{split}
\left(\bar{H}_{k} \beta e_{1}\right) =&amp;\left(\begin{array}{cccccc}
* &amp; * &amp; * &amp; * &amp; * &amp; \beta \\
* &amp; * &amp; * &amp; * &amp; * &amp; 0 \\
&amp; * &amp; * &amp; * &amp; * &amp; 0 \\
&amp; &amp; * &amp; * &amp; * &amp; 0 \\
&amp; &amp; &amp; * &amp; * &amp; 0 \\
&amp; &amp; &amp; &amp; * &amp; 0
\end{array}\right) \\
&amp;\frac{\text { Aplicando las }}{\text {rotaciones } \overrightarrow{Q_{k}}}
\\
&amp;\left(\begin{array}{cccccc}
+ &amp; + &amp; + &amp; + &amp; + &amp; \times \\
0 &amp; + &amp; + &amp; + &amp; + &amp; \times \\
&amp; 0 &amp; + &amp; + &amp; + &amp; \times \\
&amp; &amp; 0 &amp; + &amp; + &amp; \times \\
&amp; &amp; &amp; 0 &amp; + &amp; \times \\
&amp; &amp; &amp; &amp; 0 &amp; \gamma_{k}
\end{array}\right)=\left(\begin{array}{cc}
R_{n} &amp; g_{n} \\
0 &amp; \gamma_{n}
\end{array}\right)
\end{split}\]</span></p>
<p>Sujeto a que <span class="math inline">\(g_n^T \in R^n\)</span>. El algoritmo finaliza cuando <span class="math inline">\(\gamma_k &lt; \varepsilon\)</span>. <span class="citation" data-cites="Bi3"></span> <span class="citation" data-cites="Bi4"></span> <span class="citation" data-cites="Bi5"></span>. Finalmente, el algoritmo después de hallar los vector ortonormales con Arnoldi quedaría:<br />
<strong><span style="color: NavyBlue">Método Estándar (GMRES)</span></strong><br />
<br />
ENTRADA : <span class="math inline">\(x_0,bA\)</span><br />
SALIDA la solución aproximada <span class="math inline">\(x_1,...,x_n\)</span> y el residuo <span class="math inline">\(r_1,..,r_n\)</span> o un mensaje de que se excedió el numero de iteraciones.<br />
<em>Paso 1</em> Determine <span class="math inline">\(\mathbf{r_0=b-A x_0}\)</span> y <span class="math inline">\(v_1 = \mathbf{r_0/||r_0||_2}\)</span></p>
<p><em>Paso 2</em> Para <span class="math inline">\(i=1,2,...,m\)</span> haga los pasos 3, 4, 5, y 6</p>
<p><em>Paso 3</em> <span class="math inline">\(v_{j+1} = \prod_j A v_j / ||A v_j||\)</span></p>
<p><em>Paso 4</em> <span class="math inline">\(\begin{bmatrix}R_j \\ 0 \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1\bar H_j)\)</span></p>
<p><em>Paso 5</em> <span class="math inline">\(\begin{bmatrix}g_j \\ \gamma_j \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1(\beta e_1))\)</span></p>
<p><em>Paso 6</em> Si <span class="math inline">\(\gamma_j &lt; T O L\)</span>, break y haga el paso 7 y 8</p>
<p><em>Paso 7</em> <span class="math inline">\(y_k = R_k^{-1}g_k\)</span> y <span class="math inline">\(x_k = x_0 + v_k y_k\)</span></p>
<p><em>Paso 8</em> Si <span class="math inline">\(x_k\)</span> no es solución, entonces determine <span class="math inline">\(x_0=x_k\)</span> y retome el paso 1.<br />
Observe que el paso 1, y 3 esta mejor explicado en el primer pseudo-código, por lo que este ultimo se centra en la ultima parte del algoritmo que es encontrar <span class="math inline">\(y_m\)</span> y finalmente obtener <span class="math inline">\(x\)</span>.<br />
Además, podemos obtener las siguientes afirmaciones para los GMRES estándar [2, teorema 4]:</p>
<ul>
<li><p>El rango de <span class="math inline">\(A V_k\)</span> es igual al rango de <span class="math inline">\(R_k\)</span>, en particular, si <span class="math inline">\(r_{k,k}=0\)</span>, entonces A es singular.</p></li>
<li><p>El vector <span class="math inline">\(y_k\)</span> que minimiza <span class="math inline">\(\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\)</span> y <span class="math inline">\(x_n=x_0+V_ny_n\)</span>, esta dado por <span class="math inline">\(y_k=R_k^{-1}g_k\)</span></p></li>
<li><p>la normal del residual en el paso <span class="math inline">\(k\)</span> satisface <span class="math inline">\(||b - A x_k|| = |\gamma_k|\)</span></p></li>
</ul>
<h1 id="complejidades">Complejidades</h1>
<p>La complejidad computacional de GMRES <span class="math inline">\(O(n^2+n k+k^2)\)</span> y de memoria es <span class="math inline">\(O(n^2+n k+k^2)\)</span> donde el <span class="math inline">\(k\)</span> el numero de pasos de ortogonalización de. Mientras lo que se conoce del algoritmo Biconjugado es <span class="math inline">\(O(n)\)</span> en cada iteración.</p>
</body>