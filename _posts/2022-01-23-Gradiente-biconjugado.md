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
# Métodos del Gradiente Biconjugado y Mínimos Cuadrados Generalizados GMRES {#métodos-del-gradiente-biconjugado-y-mínimos-cuadrados-generalizados-gmres .title}

Miguel Gutierrez

David Felipe Martinez
:::

# Gradiente Conjugado

El método del gradiente biconjugado se basa en el método del gradiente
conjugado donde se quiere solucionar el sistema [\\\[Ax = b\\\]]{.math
.display}. Este método calcula un pseudo-gradiente [\\(\\hat
r_k\\)]{.math .inline} y una pseudo-dirección de descenso [\\(\\hat
p_k\\)]{.math .inline}. Estos pseudo gradientes serán ortogonales a los
gradientes [\\(r_k\\)]{.math .inline} y [\\(p_k\\)]{.math .inline}. Este
método a diferencia del gradiente conjugado no garantiza convergencia en
[\\(m\\)]{.math .inline} pasos.\
El pre-condicionamiento mejora en general la convergencia de los
métodos, esto pues transforma la matriz de coeficientes en una con mejor
espectro.

**[Método de gradiente biconjugado
precondicionado]{style="color: NavyBlue"}**\
ENTRADA : el numero de ecuaciones y valores desconocidos
[\\(n;\\)]{.math .inline} las entradas [\\(a\_{ij}, 1\\leq i,j\\leq
n\\)]{.math .inline} de la matriz [\\(A\\)]{.math .inline}; las entradas
[\\(b_j,1\\leq j \\leq n\\)]{.math .inline} del vector
[\\(\\mathbf{b}\\)]{.math .inline}; las entradas [\\(\\gamma\_{ij},
1\\leq i,j\\leq n\\)]{.math .inline} de la matriz precondicionada
[\\(c\^{-1}\\)]{.math .inline}, las entradas [\\(x_i,1\\leq i\\leq
n\\)]{.math .inline} de la aproximación inicial
[\\(\\mathbf{x=x\^{(0)}}\\)]{.math .inline}, el numero máximo de
iteraciones [\\(N\\)]{.math .inline}, la tolerancia [\\(TOL\\)]{.math
.inline}.\
SALIDA la solución aproximada [\\(x_1,\...,x_n\\)]{.math .inline} y el
residuo [\\(r_1,..,r_n\\)]{.math .inline} o un mensaje de que se excedió
el numero de iteraciones.\
*Paso 1* Determine [\\(\\mathbf{r=b-Ax;}\\)]{.math .inline} (*Calcule*
[\\(\\mathbf{r\^{(0)}}\\)]{.math .inline})

[\\\[\\begin{split} \\mathbf{\\hat r} &\\mathbf{= r} \\text{;
(\\textit{Calcule} \$\\mathbf{\\hat r\^{(0)}}\$)} \\\\ \\mathbf{q}
&\\mathbf{=C\^{-1}r} \\text{; (\\textit{Nota}
\$\\mathbf{q=q\^{(0)}}\$)}\\\\ \\mathbf{\\hat q}
&\\mathbf{=(C\^T)\^{-1}\\hat r}\\text{; (\\textit{Nota} \$\\mathbf{\\hat
q= \\hat q\^{(0)}}\$)} \\\\ \\end{split}\\\]]{.math .display}

*Paso 2* Determine [\\(k=1\\)]{.math .inline}

*Paso 3* Mientras ([\\(k \\leq N\\)]{.math .inline}) haga los pasos 4-7

*Paso 4* Si [\\(\|\|r_k\|\| \\leq TOL\\)]{.math .inline}, entonces,
retornar.

*Paso 5* Determine [\\\[\\begin{split} \\alpha &= \\mathbf{\\frac{\\hat
r \\cdot q}{\\hat p\\cdot A \\cdot p}} \\text{; (\\textit{Nota}
\$\\mathbf{\\hat \\alpha= \\frac{ \\langle \\hat r\^{(k)},
q\^{(k)}\\rangle} {\\langle \\hat p\^{(k)},A p\^{(k)} \\rangle}}\$)}\\\\
\\mathbf{r} &\\mathbf{=r-\\alpha q} \\text{; (\\textit{Nota}
\$\\mathbf{r=r\^{(k)}}\$)}\\\\ \\hat \\alpha &= \\mathbf{\\frac{\\hat r
\\cdot \\hat q}{\\hat p\\cdot A\^T \\cdot \\hat p}} \\text{;
(\\textit{Nota} \$\\mathbf{\\hat \\alpha= \\frac{ \\langle \\hat
r\^{(k)}, \\hat q\^{(k)}\\rangle} {\\langle \\hat p\^{(k)},A\^T \\hat
p\^{(k)} \\rangle}}\$)}\\\\ \\mathbf{\\hat r} & \\mathbf{=\\hat r-\\hat
\\alpha \\hat q} \\\\ \\mathbf{x\_{k+1}} &\\mathbf{= x_k + \\alpha p}
\\\\ \\end{split}\\\]]{.math .display}

\
*Paso 6* Determine [\\\[\\begin{split} \\mathbf{q} & \\mathbf{=
C\^{-1}r} \\text{; (\\textit{Nota} \$\\mathbf{q=q\^{(k)}}\$)}\\\\
\\mathbf{\\hat q} &\\mathbf{= (C\^T)\^{-1}\\hat r} \\text{;
(\\textit{Nota} \$\\mathbf{\\hat q=\\hat q\^{(k)}}\$)} \\\\
\\mathbf{\\beta} &\\mathbf{= \\frac{\\hat r \\cdot q}{\\hat r \\cdot
q\_{k-1}}} \\text{; (\\textit{Nota} \$\\mathbf{\\beta= \\frac{ \\langle
\\hat r\^{(k)}, q\^{(k+1)}\\rangle} {\\langle \\hat r\^{(k)}, q\^{(k)}
\\rangle}}\$)} \\end{split}\\\]]{.math .display}

*Paso 7* Determine [\\\[\\begin{split} \\mathbf{p} &\\mathbf{=
q\_{k+1}+\\beta p} \\text{; (\\textit{Nota}
\$\\mathbf{p=p\^{(k)}}\$)}\\\\ \\mathbf{\\hat p} &\\mathbf{= \\hat
q\_{k+1}+\\beta \\hat p} \\text{; (\\textit{Nota} \$\\mathbf{\\hat p=
\\hat p\^{(k)}}\$)} \\\\ k &= k+1 \\end{split}\\\]]{.math .display}

*Paso 8* Si [\\((k>n)\\)]{.math .inline} entonces SALIDA. PARE.\
La secuencia de vectores satisfacen la condición de biortogonalidad

[\\\[\\mathbf{\\hat r_i \\cdot r_j = r_i \\cdot \\hat r_j}= 0, \\quad
j\<i\\\]]{.math .display}

y la condición de biconjugancia

[\\\[\\mathbf{\\hat p_i \\cdot A \\cdot p_j = p_i \\cdot A\^T \\cdot
\\hat p_j}= 0, \\quad j\<i\\\]]{.math .display}

también se cumple ortogonalidad mutua

[\\\[\\mathbf{\\hat r_i \\cdot p_j = r_i \\cdot \\hat p_j}= 0, \\quad
j\<i\\\]]{.math .display}

Las demostraciones de estas propiedades resultan del proceso de
inducción. Note que cuando [\\(A\\)]{.math .inline} es simétrica y se
escoge [\\(\\hat r = r\\)]{.math .inline}, luego [\\(\\hat r_k =
r_k\\)]{.math .inline} y [\\(\\hat p_k = p_k\\)]{.math .inline} para
todo [\\(k\\)]{.math .inline} es el caso del algoritmo del gradiente
precondicionado. La matriz condicionadora mejora la convergencia del
algoritmo []{.citation cites="Bi1"}. Note que el gradiente conjugado es
el doble de costoso que el método del gradiente.

# Residuo Mínimo Generalizado (GMRES)

El método del Residuo Mínimo Generalizado (GMRES) es un método iterativo
utilizado para encontrar la soluciona a un sistema lineal no simétrico.
El método reside en reconstruir una base ortonormal del espacio de
Krylov y por lo tanto es vulnerable a sesgos de errores computacionales.
Este método reduce las dimensiones y por lo tanto resulta eficiente
comparado con otros métodos.

El espacio de Krylov de dimensión [\\(k\\)]{.math .inline} se define
como el subespacio lineal mas pequeño que contiene el conjunto por las
imágenes de [\\(b\\)]{.math .inline} bajo las primeras [\\(k\\)]{.math
.inline} potencias de A.

[\\\[\\begin{split} \\kappa_n &= span\\{b,Ab,\...,A\^{n-1}b\\}\\\\
&=\\{q_1,q2,..,q_n\\} \\subset \\mathbf{R}\^{m} \\end{split}\\\]]{.math
.display}

Donde cada [\\(q_i\\)]{.math .inline} son vectores ortonormales. Ahora
definimos la matriz de Krylov:

[\\\[K_n = \\begin{bmatrix} \\vdots & \\vdots &\\vdots &\\vdots\\\\ b &
Ab & A\^2b & A\^{n-1}b \\\\ \\vdots & \\vdots &\\vdots &\\vdots
\\end{bmatrix} = V_n R_n\\\]]{.math .display} Donde [\\(V_n\\)]{.math
.inline} y [\\(R_n\\)]{.math .inline} conforman las matrices de una
descomposición [\\(Q R\\)]{.math .inline} respectivamente.\
Cualquier vector [\\(x\\in K_n\\)]{.math .inline} puede ser escrito como
[\\(x = x_0 + V_n y_n\\)]{.math .inline} donde [\\(y\\)]{.math .inline}
es un [\\(n\\)]{.math .inline}-vector y [\\(V_n\\)]{.math .inline} es
una base ortonormal del subespacio de Krylov [\\(K_n\\)]{.math .inline}.
Definimos la función de perdida que queremos minimizar.

[\\\[\\begin{split} J(y) = \|\|b-Ax\|\|\_2 &=\|\|b-A(x_0+V_n
y)\|\|\_2\\\\ &= \|\|r\^0-AV_n y\|\|\_2 \\\\ \\end{split}\\\]]{.math
.display}

Con la relación de arnoldi [\\\[AV_n = V\_{n+1}\\bar H_k\\\]]{.math
.display}

obtenemos

[\\\[\\begin{split} \\left\\\|r\_{0}-A V\_{n}
y\\right\\\|\_{2}&=\\left\\\|r\_{0}-V\_{n+1} \\bar{H}\_{n}
y\\right\\\|\_{2}\\\\ &=\\left\\\|V\_{n+1}\\left(\\beta
e\_{1}-\\bar{H}\_{n} y\\right)\\right\\\|\_{2}\\\\ &=\\left\\\|\\beta
e\_{1}-\\bar{H}\_{n} y\\right\\\|\_{2} \\end{split}\\\]]{.math .display}

Donde [\\(r_0 = \\beta v_1\\)]{.math .inline}, [\\(e_1\\)]{.math
.inline} es la primera columna de la matriz identidad
[\\((k+1)\\times(k+1)\\)]{.math .inline} y [\\(\\bar H\\)]{.math
.inline} es una matriz superior de Hessenberg

[\\\[\\bar H\_{n}=\\left\[\\begin{array}{cccc} h\_{11} & h\_{12} &
\\cdots & h\_{1 n} \\\\ h\_{21} & h\_{22} & \\cdots & h\_{2 n} \\\\ 0 &
h\_{32} & \\ddots & h\_{3 n} \\\\ \\vdots & \\vdots & & \\vdots \\\\ 0 &
0 & \\cdots & h\_{n n} \\end{array}\\right\]\\\]]{.math .display}

Luego el problema se convierte en minimizar en cada iteración

[\\\[\\begin{split} x\_{n}&=x\_{0}+V\_{n} y\_{n} \\\\ y\_{n}&=\\arg
\\min \_{y}\\left\\\|\\left(\\beta e\_{1}-\\bar{H}\_{n}
y\\right)\\right\\\|\_{2} \\end{split}\\\]]{.math .display}

Para encontrar [\\(V_n\\)]{.math .inline} el algoritmo de Arnoldi nos
será de utilidad.\
**[Método de residuo Mínimo Generalizado
(GMRES)]{style="color: NavyBlue"}**\
ENTRADA : [\\(x_0,b,A\\)]{.math .inline}\
SALIDA la solución aproximada [\\(x_1,\...,x_n\\)]{.math .inline} y el
residuo [\\(r_1,..,r_n\\)]{.math .inline} o un mensaje de que se excedió
el numero de iteraciones.\
*Paso 1* Determine [\\(\\mathbf{r_0=b-Ax_0}\\)]{.math .inline} y [\\(v_1
= \\mathbf{r_0/\|\|r_0\|\|\_2}\\)]{.math .inline}

*Paso 2* Defina una matriz superior de Hesenberg [\\(H_n=0\\)]{.math
.inline} de tamaño [\\((m+1)\\times m\\)]{.math .inline}.

*Paso 3* Para [\\(j=1,2,\...,m\\)]{.math .inline} haga los pasos 4, 5, 6
y 7

*Paso 4* Determine [\\(w_j = Av_j\\)]{.math .inline}

*Paso 5* Para [\\(i=1,2,\...,m\\)]{.math .inline} determine
[\\\[\\begin{split} h\_{ij} &= \\langle w_j,vi \\rangle \\\\ w_j &= w_j
- h\_{ij}v_i \\end{split}\\\]]{.math .display}

*Paso 6* Determine [\\(h\_{j+1j} = \|\|w_j\|\|\_2\\)]{.math .inline}.

*Paso 7* Si [\\(h\_{j+1,j}=0\\)]{.math .inline} entonces
[\\(m=j\\)]{.math .inline} y romper iteración (break).

*Paso 8* Determine [\\(v\_{j+1} = w_j/h\_{j+1,j}\\)]{.math .inline}.

*Paso 9* Computar minimizador de [\\(\\left\\\|\\left(\\beta
e\_{1}-\\bar{H}\_{n} y\\right)\\right\\\|\_{2}\\)]{.math .inline} y
[\\(x_n=x_0+V_ny_n\\)]{.math .inline}.\
Note que en los pasos [\\(1-8\\)]{.math .inline} se esta encontrando una
base ortonormal del espacio de Krylov. De aquí en adelante, la pregunta
¿Como encontramos [\\(y_m\\)]{.math .inline} ? Es donde surgen
diferentes métodos. Uno de los mas simple y poderosos fue introducido
por Saad Y Schultz.

## GMRES Estandar

Para solucionar el problema de mínimos cuadrados

[\\\[y\_{n} =\\arg \\min \_{y}\\left\\\|\\left(\\beta
e\_{1}-\\bar{H}\_{n} y\\right)\\right\\\|\_{2}\\\]]{.math .display}

Se sugirió utilizar unas rotaciones dadas para convertir la anterior
ecuación en un sistema superior triangular para descomponer
[\\(H_n\\)]{.math .inline} en [\\(Q_n \\hat R_k\\)]{.math .inline} Donde
[\\(Q_n\\)]{.math .inline} es una matriz ortonormal obtenida del
producto de [\\(k\\)]{.math .inline} rotaciones y [\\(\\hat R_k
=\\begin{bmatrix}R_k \\\\ 0\...0 \\end{bmatrix}\\)]{.math .inline}. La
matriz de rotación es una matriz identidad de orden [\\(n+1\\)]{.math
.inline} donde solo cuatro componentes son remplazadas por escalares
[\\(c_i,s_i\\)]{.math .inline} de la siguiente manera

[\\\[J_i = \\begin{bmatrix} 1 & & & 0 \\\\ &c_i & s_i & \\\\ & -s_i &c_i
& \\\\ 0 & & & 1 \\\\ \\end{bmatrix} \\quad \\text{donde } c_i\^2+s_i\^2
= 1,i= 1,\...,n\\\]]{.math .display}

La matriz es construida de tal manera que

[\\\[\\begin{bmatrix} c_i & s_i \\\\ -s_i &c_i \\\\ \\end{bmatrix}
\\begin{bmatrix} h\_{i,i} \\\\ h\_{i+1,i} \\\\ \\end{bmatrix} =
\\begin{bmatrix} \* \\\\ 0 \\\\ \\end{bmatrix}\\\]]{.math .display}

El producto punto de las primeras [\\(n\\)]{.math .inline} rotaciones,
es decir, [\\(Q_n = J_nJ\_{k-1}\\hdots J_1\\)]{.math .inline} al
aplicarse sobre la matriz por izquierda de [\\(H_n\\)]{.math .inline} y
[\\(\\beta e_1\\)]{.math .inline} donde se obtiene la descomposición

[\\\[R_ny_n = g_n\\\]]{.math .display}

y la norma residual se obtiene

[\\\[\\begin{split} \\left(\\bar{H}\_{k} \\beta e\_{1}\\right)
=&\\left(\\begin{array}{cccccc} \* & \* & \* & \* & \* & \\beta \\\\ \*
& \* & \* & \* & \* & 0 \\\\ & \* & \* & \* & \* & 0 \\\\ & & \* & \* &
\* & 0 \\\\ & & & \* & \* & 0 \\\\ & & & & \* & 0 \\end{array}\\right)
\\\\ &\\frac{\\text { Aplicando las }}{\\text {rotaciones }
\\overrightarrow{Q\_{k}}} \\\\ &\\left(\\begin{array}{cccccc} + & + & +
& + & + & \\times \\\\ 0 & + & + & + & + & \\times \\\\ & 0 & + & + & +
& \\times \\\\ & & 0 & + & + & \\times \\\\ & & & 0 & + & \\times \\\\ &
& & & 0 & \\gamma\_{k} \\end{array}\\right)=\\left(\\begin{array}{cc}
R\_{n} & g\_{n} \\\\ 0 & \\gamma\_{n} \\end{array}\\right)
\\end{split}\\\]]{.math .display}

Sujeto a que [\\(g_n\^T \\in R\^n\\)]{.math .inline}. El algoritmo
finaliza cuando [\\(\\gamma_k \< \\varepsilon\\)]{.math .inline}.
[]{.citation cites="Bi3"} []{.citation cites="Bi4"} []{.citation
cites="Bi5"}. Finalmente, el algoritmo después de hallar los vector
ortonormales con Arnoldi quedaría:\
**[Método Estándar (GMRES)]{style="color: NavyBlue"}**\
\
ENTRADA : [\\(x_0,bA\\)]{.math .inline}\
SALIDA la solución aproximada [\\(x_1,\...,x_n\\)]{.math .inline} y el
residuo [\\(r_1,..,r_n\\)]{.math .inline} o un mensaje de que se excedió
el numero de iteraciones.\
*Paso 1* Determine [\\(\\mathbf{r_0=b-A x_0}\\)]{.math .inline} y
[\\(v_1 = \\mathbf{r_0/\|\|r_0\|\|\_2}\\)]{.math .inline}

*Paso 2* Para [\\(i=1,2,\...,m\\)]{.math .inline} haga los pasos 3, 4,
5, y 6

*Paso 3* [\\(v\_{j+1} = \\prod_j A v_j / \|\|A v_j\|\|\\)]{.math
.inline}

*Paso 4* [\\(\\begin{bmatrix}R_j \\\\ 0 \\\\\\end{bmatrix} =
J_j(J\_{n-1} \\cdots J_1\\bar H_j)\\)]{.math .inline}

*Paso 5* [\\(\\begin{bmatrix}g_j \\\\ \\gamma_j \\\\\\end{bmatrix} =
J_j(J\_{n-1} \\cdots J_1(\\beta e_1))\\)]{.math .inline}

*Paso 6* Si [\\(\\gamma_j \< T O L\\)]{.math .inline}, break y haga el
paso 7 y 8

*Paso 7* [\\(y_k = R_k\^{-1}g_k\\)]{.math .inline} y [\\(x_k = x_0 + v_k
y_k\\)]{.math .inline}

*Paso 8* Si [\\(x_k\\)]{.math .inline} no es solución, entonces
determine [\\(x_0=x_k\\)]{.math .inline} y retome el paso 1.\
Observe que el paso 1, y 3 esta mejor explicado en el primer
pseudo-código, por lo que este ultimo se centra en la ultima parte del
algoritmo que es encontrar [\\(y_m\\)]{.math .inline} y finalmente
obtener [\\(x\\)]{.math .inline}.\
Además, podemos obtener las siguientes afirmaciones para los GMRES
estándar \[2, teorema 4\]:

-   El rango de [\\(A V_k\\)]{.math .inline} es igual al rango de
    [\\(R_k\\)]{.math .inline}, en particular, si
    [\\(r\_{k,k}=0\\)]{.math .inline}, entonces A es singular.

-   El vector [\\(y_k\\)]{.math .inline} que minimiza
    [\\(\\left\\\|\\left(\\beta e\_{1}-\\bar{H}\_{n}
    y\\right)\\right\\\|\_{2}\\)]{.math .inline} y
    [\\(x_n=x_0+V_ny_n\\)]{.math .inline}, esta dado por
    [\\(y_k=R_k\^{-1}g_k\\)]{.math .inline}

-   la normal del residual en el paso [\\(k\\)]{.math .inline} satisface
    [\\(\|\|b - A x_k\|\| = \|\\gamma_k\|\\)]{.math .inline}

# Complejidades

La complejidad computacional de GMRES [\\(O(n\^2+n k+k\^2)\\)]{.math
.inline} y de memoria es [\\(O(n\^2+n k+k\^2)\\)]{.math .inline} donde
el [\\(k\\)]{.math .inline} el numero de pasos de ortogonalización de.
Mientras lo que se conoce del algoritmo Biconjugado es
[\\(O(n)\\)]{.math .inline} en cada iteración.
