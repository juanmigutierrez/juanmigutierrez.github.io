---
layout: post
read_time: true
show_date: true
title: "Gradiente Biconjugado"
date: 2022-01-23
img: posts/20210420/post8-rembrandt.jpg
tags: [gradiente,biconjugado]
category: summary
author: Miguel Gutierrez y David Felipe Martinez
description: "Resumen gradiente biconjugado"
---


Gradiente Conjugado
===================

El método del gradiente biconjugado se basa en el método del gradiente
conjugado donde se quiere solucionar el sistema $$Ax = b$$. Este método
calcula un pseudo-gradiente $\hat r_k$ y una pseudo-dirección de
descenso $\hat p_k$. Estos pseudo gradientes serán ortogonales a los
gradientes $r_k$ y $p_k$. Este método a diferencia del gradiente
conjugado no garantiza convergencia en $m$ pasos.\
El pre-condicionamiento mejora en general la convergencia de los
métodos, esto pues transforma la matriz de coeficientes en una con mejor
espectro.

**[Método de gradiente biconjugado
precondicionado]{style="color: NavyBlue"}**\
ENTRADA : el numero de ecuaciones y valores desconocidos $n;$ las
entradas $a_{ij}, 1\leq i,j\leq n$ de la matriz $A$; las entradas
$b_j,1\leq j \leq n$ del vector $\mathbf{b}$; las entradas
$\gamma_{ij}, 1\leq i,j\leq n$ de la matriz precondicionada $c^{-1}$,
las entradas $x_i,1\leq i\leq n$ de la aproximación inicial
$\mathbf{x=x^{(0)}}$, el numero máximo de iteraciones $N$, la tolerancia
$TOL$.\
SALIDA la solución aproximada $x_1,...,x_n$ y el residuo $r_1,..,r_n$ o
un mensaje de que se excedió el numero de iteraciones.\
*Paso 1* Determine $\mathbf{r=b-Ax;}$ (*Calcule* $\mathbf{r^{(0)}}$)

$$\begin{split}
   \mathbf{\hat r} &\mathbf{= r} \text{; (\textit{Calcule} $\mathbf{\hat r^{(0)}}$)} \\
   \mathbf{q} &\mathbf{=C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(0)}}$)}\\
   \mathbf{\hat q} &\mathbf{=(C^T)^{-1}\hat r}\text{; (\textit{Nota} $\mathbf{\hat q= \hat q^{(0)}}$)} \\
\end{split}$$

*Paso 2* Determine $k=1$

*Paso 3* Mientras ($k \leq N$) haga los pasos 4-7

*Paso 4* Si $||r_k|| \leq TOL$, entonces, retornar.

*Paso 5* Determine $$\begin{split}
   \alpha &= \mathbf{\frac{\hat r \cdot q}{\hat p\cdot A \cdot p}}  \text{; (\textit{Nota} $\mathbf{\hat \alpha=  \frac{ \langle \hat r^{(k)},  q^{(k)}\rangle} {\langle \hat p^{(k)},A p^{(k)} \rangle}}$)}\\
   \mathbf{r} &\mathbf{=r-\alpha q} \text{; (\textit{Nota} $\mathbf{r=r^{(k)}}$)}\\
    \hat \alpha &= \mathbf{\frac{\hat r \cdot \hat q}{\hat p\cdot A^T \cdot \hat p}} \text{; (\textit{Nota} $\mathbf{\hat \alpha=  \frac{ \langle \hat r^{(k)}, \hat q^{(k)}\rangle} {\langle \hat p^{(k)},A^T \hat p^{(k)} \rangle}}$)}\\
   \mathbf{\hat r} & \mathbf{=\hat r-\hat \alpha \hat q} \\
   \mathbf{x_{k+1}} &\mathbf{= x_k + \alpha p} \\
\end{split}$$

\
*Paso 6* Determine $$\begin{split}
   \mathbf{q} & \mathbf{= C^{-1}r} \text{; (\textit{Nota} $\mathbf{q=q^{(k)}}$)}\\
   \mathbf{\hat q} &\mathbf{= (C^T)^{-1}\hat r} \text{; (\textit{Nota} $\mathbf{\hat q=\hat q^{(k)}}$)} \\
   \mathbf{\beta} &\mathbf{= \frac{\hat r \cdot q}{\hat r \cdot q_{k-1}}} \text{; (\textit{Nota} $\mathbf{\beta=  \frac{ \langle \hat r^{(k)}, q^{(k+1)}\rangle} {\langle \hat r^{(k)},   q^{(k)} \rangle}}$)}
\end{split}$$

*Paso 7* Determine $$\begin{split}
   \mathbf{p} &\mathbf{= q_{k+1}+\beta p} \text{; (\textit{Nota} $\mathbf{p=p^{(k)}}$)}\\
   \mathbf{\hat p} &\mathbf{= \hat q_{k+1}+\beta \hat p} \text{; (\textit{Nota} $\mathbf{\hat p= \hat p^{(k)}}$)} \\
   k &= k+1
\end{split}$$

*Paso 8* Si $(k>n)$ entonces SALIDA. PARE.\
La secuencia de vectores satisfacen la condición de biortogonalidad

$$\mathbf{\hat r_i \cdot r_j = r_i \cdot \hat r_j}= 0, \quad j<i$$

y la condición de biconjugancia

$$\mathbf{\hat p_i \cdot A \cdot p_j = p_i \cdot A^T \cdot \hat p_j}= 0, \quad j<i$$

también se cumple ortogonalidad mutua

$$\mathbf{\hat r_i \cdot p_j = r_i \cdot \hat p_j}= 0, \quad j<i$$

Las demostraciones de estas propiedades resultan del proceso de
inducción. Note que cuando $A$ es simétrica y se escoge $\hat r = r$,
luego $\hat r_k = r_k$ y $\hat p_k = p_k$ para todo $k$ es el caso del
algoritmo del gradiente precondicionado. La matriz condicionadora mejora
la convergencia del algoritmo [@Bi1]. Note que el gradiente conjugado es
el doble de costoso que el método del gradiente.

Residuo Mínimo Generalizado (GMRES)
===================================

El método del Residuo Mínimo Generalizado (GMRES) es un método iterativo
utilizado para encontrar la soluciona a un sistema lineal no simétrico.
El método reside en reconstruir una base ortonormal del espacio de
Krylov y por lo tanto es vulnerable a sesgos de errores computacionales.
Este método reduce las dimensiones y por lo tanto resulta eficiente
comparado con otros métodos.

El espacio de Krylov de dimensión $k$ se define como el subespacio
lineal mas pequeño que contiene el conjunto por las imágenes de $b$ bajo
las primeras $k$ potencias de A.

$$\begin{split}
    \kappa_n &= span\{b,Ab,...,A^{n-1}b\}\\
    &=\{q_1,q2,..,q_n\}  \subset \mathbf{R}^{m}
\end{split}$$

Donde cada $q_i$ son vectores ortonormales. Ahora definimos la matriz de
Krylov:

$$K_n =
\begin{bmatrix}
\vdots & \vdots &\vdots &\vdots\\
b & Ab & A^2b & A^{n-1}b \\
\vdots & \vdots &\vdots &\vdots
\end{bmatrix}
= V_n R_n$$ Donde $V_n$ y $R_n$ conforman las matrices de una
descomposición $Q R$ respectivamente.\
Cualquier vector $x\in K_n$ puede ser escrito como $x = x_0 + V_n y_n$
donde $y$ es un $n$-vector y $V_n$ es una base ortonormal del subespacio
de Krylov $K_n$. Definimos la función de perdida que queremos minimizar.

$$\begin{split}
    J(y) = ||b-Ax||_2 &=||b-A(x_0+V_n y)||_2\\
    &= ||r^0-AV_n y||_2 \\
\end{split}$$

Con la relación de arnoldi $$AV_n = V_{n+1}\bar H_k$$

obtenemos

$$\begin{split}
\left\|r_{0}-A V_{n} y\right\|_{2}&=\left\|r_{0}-V_{n+1} \bar{H}_{n} y\right\|_{2}\\
&=\left\|V_{n+1}\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}\\
&=\left\|\beta e_{1}-\bar{H}_{n} y\right\|_{2}
\end{split}$$

Donde $r_0 = \beta v_1$, $e_1$ es la primera columna de la matriz
identidad $(k+1)\times(k+1)$ y $\bar H$ es una matriz superior de
Hessenberg

$$\bar H_{n}=\left[\begin{array}{cccc}
h_{11} & h_{12} & \cdots & h_{1 n} \\
h_{21} & h_{22} & \cdots & h_{2 n} \\
0 & h_{32} & \ddots & h_{3 n} \\
\vdots & \vdots & & \vdots \\
0 & 0 & \cdots & h_{n n}
\end{array}\right]$$

Luego el problema se convierte en minimizar en cada iteración

$$\begin{split}
x_{n}&=x_{0}+V_{n} y_{n} \\
y_{n}&=\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}
\end{split}$$

Para encontrar $V_n$ el algoritmo de Arnoldi nos será de utilidad.\
**[Método de residuo Mínimo Generalizado
(GMRES)]{style="color: NavyBlue"}**\
ENTRADA : $x_0,b,A$\
SALIDA la solución aproximada $x_1,...,x_n$ y el residuo $r_1,..,r_n$ o
un mensaje de que se excedió el numero de iteraciones.\
*Paso 1* Determine $\mathbf{r_0=b-Ax_0}$ y
$v_1 = \mathbf{r_0/||r_0||_2}$

*Paso 2* Defina una matriz superior de Hesenberg $H_n=0$ de tamaño
$(m+1)\times m$.

*Paso 3* Para $j=1,2,...,m$ haga los pasos 4, 5, 6 y 7

*Paso 4* Determine $w_j = Av_j$

*Paso 5* Para $i=1,2,...,m$ determine $$\begin{split}
    h_{ij} &= \langle w_j,vi \rangle \\
    w_j &= w_j - h_{ij}v_i
\end{split}$$

*Paso 6* Determine $h_{j+1j} = ||w_j||_2$.

*Paso 7* Si $h_{j+1,j}=0$ entonces $m=j$ y romper iteración (break).

*Paso 8* Determine $v_{j+1} = w_j/h_{j+1,j}$.

*Paso 9* Computar minimizador de
$\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}$ y
$x_n=x_0+V_ny_n$.\
Note que en los pasos $1-8$ se esta encontrando una base ortonormal del
espacio de Krylov. De aquí en adelante, la pregunta ¿Como encontramos
$y_m$ ? Es donde surgen diferentes métodos. Uno de los mas simple y
poderosos fue introducido por Saad Y Schultz.

GMRES Estandar
--------------

Para solucionar el problema de mínimos cuadrados

$$y_{n} =\arg \min _{y}\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}$$

Se sugirió utilizar unas rotaciones dadas para convertir la anterior
ecuación en un sistema superior triangular para descomponer $H_n$ en
$Q_n \hat R_k$ Donde $Q_n$ es una matriz ortonormal obtenida del
producto de $k$ rotaciones y $\hat R_k =\begin{bmatrix}R_k \\ 0...0
\end{bmatrix}$. La matriz de rotación es una matriz identidad de orden
$n+1$ donde solo cuatro componentes son remplazadas por escalares
$c_i,s_i$ de la siguiente manera

$$J_i =
\begin{bmatrix}
1 & & & 0 \\
&c_i & s_i & \\ 
& -s_i &c_i & \\
0 & & & 1 \\
\end{bmatrix}
\quad \text{donde } c_i^2+s_i^2 = 1,i= 1,...,n$$

La matriz es construida de tal manera que

$$\begin{bmatrix}
c_i & s_i  \\ 
-s_i &c_i \\
\end{bmatrix}
\begin{bmatrix}
h_{i,i}  \\ 
h_{i+1,i} \\
\end{bmatrix}
= 
\begin{bmatrix}
* \\ 
0 \\
\end{bmatrix}$$

El producto punto de las primeras $n$ rotaciones, es decir,
$Q_n = J_nJ_{k-1}\hdots J_1$ al aplicarse sobre la matriz por izquierda
de $H_n$ y $\beta e_1$ donde se obtiene la descomposición

$$R_ny_n = g_n$$

y la norma residual se obtiene

$$\begin{split}
\left(\bar{H}_{k} \beta e_{1}\right) =&\left(\begin{array}{cccccc}
* & * & * & * & * & \beta \\
* & * & * & * & * & 0 \\
& * & * & * & * & 0 \\
& & * & * & * & 0 \\
& & & * & * & 0 \\
& & & & * & 0
\end{array}\right) \\
&\frac{\text { Aplicando las }}{\text {rotaciones } \overrightarrow{Q_{k}}}
\\
&\left(\begin{array}{cccccc}
+ & + & + & + & + & \times \\
0 & + & + & + & + & \times \\
& 0 & + & + & + & \times \\
& & 0 & + & + & \times \\
& & & 0 & + & \times \\
& & & & 0 & \gamma_{k}
\end{array}\right)=\left(\begin{array}{cc}
R_{n} & g_{n} \\
0 & \gamma_{n}
\end{array}\right)
\end{split}$$

Sujeto a que $g_n^T \in R^n$. El algoritmo finaliza cuando
$\gamma_k < \varepsilon$. [@Bi3] [@Bi4] [@Bi5]. Finalmente, el algoritmo
después de hallar los vector ortonormales con Arnoldi quedaría:\
**[Método Estándar (GMRES)]{style="color: NavyBlue"}**\
\
ENTRADA : $x_0,bA$\
SALIDA la solución aproximada $x_1,...,x_n$ y el residuo $r_1,..,r_n$ o
un mensaje de que se excedió el numero de iteraciones.\
*Paso 1* Determine $\mathbf{r_0=b-A x_0}$ y
$v_1 = \mathbf{r_0/||r_0||_2}$

*Paso 2* Para $i=1,2,...,m$ haga los pasos 3, 4, 5, y 6

*Paso 3* $v_{j+1} = \prod_j A v_j / ||A v_j||$

*Paso 4*
$\begin{bmatrix}R_j \\ 0 \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1\bar H_j)$

*Paso 5*
$\begin{bmatrix}g_j \\ \gamma_j \\\end{bmatrix} = J_j(J_{n-1} \cdots J_1(\beta e_1))$

*Paso 6* Si $\gamma_j < T O L$, break y haga el paso 7 y 8

*Paso 7* $y_k = R_k^{-1}g_k$ y $x_k = x_0 + v_k y_k$

*Paso 8* Si $x_k$ no es solución, entonces determine $x_0=x_k$ y retome
el paso 1.\
Observe que el paso 1, y 3 esta mejor explicado en el primer
pseudo-código, por lo que este ultimo se centra en la ultima parte del
algoritmo que es encontrar $y_m$ y finalmente obtener $x$.\
Además, podemos obtener las siguientes afirmaciones para los GMRES
estándar \[2, teorema 4\]:

-   El rango de $A V_k$ es igual al rango de $R_k$, en particular, si
    $r_{k,k}=0$, entonces A es singular.

-   El vector $y_k$ que minimiza
    $\left\|\left(\beta e_{1}-\bar{H}_{n} y\right)\right\|_{2}$ y
    $x_n=x_0+V_ny_n$, esta dado por $y_k=R_k^{-1}g_k$

-   la normal del residual en el paso $k$ satisface
    $||b - A x_k|| = |\gamma_k|$

Complejidades
=============

La complejidad computacional de GMRES $O(n^2+n k+k^2)$ y de memoria es
$O(n^2+n k+k^2)$ donde el $k$ el numero de pasos de ortogonalización de.
Mientras lo que se conoce del algoritmo Biconjugado es $O(n)$ en cada
iteración.
