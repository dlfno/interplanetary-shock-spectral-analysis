# Fundamento científico — análisis espectral de un choque no colisional

> Notas complementarias del repositorio **Análisis espectral de un choque
> interplanetario**. Amplían el resumen condensado del `README.md`.
> La versión en inglés está en [`theory_en.md`](theory_en.md).
>
> *(Las ecuaciones están en LaTeX compatible con GitHub y se renderizan en github.com.)*

---

## (a) Choques no colisionales e interplanetarios

Un choque es una transición irreversible que **frena** un flujo
super-magnetosónico a velocidad sub-magnetosónica mientras **comprime y
calienta** el plasma. En un gas ordinario la disipación la aportan las colisiones
de Coulomb y el frente tiene un grosor de unos pocos recorridos libres medios. En
el viento solar el recorrido libre medio de Coulomb es del orden de **1 UA**
($\sim 10^8$ km), y sin embargo los frentes de choque observados solo tienen
$\sim 10^2\text{-}10^3$ km de grosor. La disipación debe ser, por tanto,
**no colisional** — proporcionada por procesos electromagnéticos colectivos:
interacciones onda–partícula, micro-inestabilidades y reflexión de iones (una
"resistividad/viscosidad anómala") que termalizan el flujo a escalas cinéticas.
Las longitudes relevantes son la longitud inercial iónica $l_i = c/\omega_{pi}$
($\approx 228/\sqrt{n\,[\mathrm{cm^{-3}}]}$ km, $\approx 130$ km con
$n = 3.1\ \mathrm{cm^{-3}}$) y el radio de giro iónico $r_{ci} = v_{th}/\Omega_{ci}$.

La MHD ideal admite tres modos de onda: los modos compresivos **rápido** y
**lento** magnetosónicos y el modo intermedio no compresivo de **Alfvén**
($v_A = B/\sqrt{\mu_0\rho}$). Las velocidades de fase magnetosónicas son

```math
v_{f,s}^2 = \tfrac{1}{2}\Big[(v_A^2 + c_s^2) \pm \sqrt{(v_A^2 + c_s^2)^2 - 4\,v_A^2 c_s^2 \cos^2\theta}\,\Big],
\qquad c_s = \sqrt{\gamma p/\rho}.
```

Los choques interplanetarios (IP) son mayoritariamente choques **magnetosónicos
rápidos**. Un choque **rápido directo (FF, *fast-forward*)** se propaga en
sentido anti-solar más rápido que el viento corriente arriba (el tipo más común a
1 UA, generalmente impulsado por ICMEs o regiones de interacción de flujos);
**B y la densidad aumentan a través del frente** — justo nuestro evento del
2022-08-19. Las relaciones de **Rankine–Hugoniot** conservan masa, momento,
energía, $\nabla\cdot\mathbf{B}=0$ y la condición de congelamiento del flujo a
través del frente. En el límite gasdinámico la compresión es
$r = \rho_d/\rho_u = (\gamma+1)M^2/[(\gamma-1)M^2 + 2]$, con máximo de choque
fuerte $r = (\gamma+1)/(\gamma-1) = 4$ para $\gamma = 5/3$; la compresión
magnética $B_d/B_u$ también es $\le 4$. Los choques reales a 1 UA tienen
$r \approx 1.5\text{-}3$ — nuestro evento tiene una compresión **moderada**
$r \approx 2.4$ (ajuste CfA), consistente con el salto observado de $|B|$ de
$\approx 5.5$ a $\approx 11$ nT.

El **ángulo normal del choque** $\theta_{Bn}$ (entre la $\mathbf{B}$ corriente
arriba y la normal $\mathbf{n}$ del choque) divide los choques en
**casi-perpendiculares** ($\theta_{Bn} > 45^\circ$; delgados, laminares, con
pie–rampa–sobrepaso) y **casi-paralelos** ($\theta_{Bn} < 45^\circ$; extendidos,
turbulentos, con iones reflejados que fluyen lejos corriente arriba formando un
*antechoque* (foreshock) y SLAMS). Con $\theta_{Bn} \approx 31^\circ$ nuestro
choque es **casi-paralelo**. Un choque requiere número de Mach rápido
$M_f = v_{un}/v_f > 1$; el **primer número de Mach crítico**
($1 \le M_c \le 2.76$) separa los choques subcríticos (resistivos) de los
supercríticos (dominados por reflexión de iones). Con $M_f \approx 1.81$ es un
choque FF casi-paralelo débil-a-moderado, cercano al crítico.

## (b) Actividad de ondas corriente arriba vs corriente abajo

La geometría fija la física de las ondas. Los choques casi-paralelos permiten que
los iones escapen hacia el Sol a lo largo del campo, produciendo un **antechoque
iónico extendido** (haces alineados al campo → intermedios → en giro → difusos)
lleno de ondas. El modo dominante es la **onda ULF de "30 segundos"**:
casi-monocromática, de gran amplitud ($\delta B/B \sim 0.1\text{-}1$),
magnetosónica rápida, a una frecuencia en el marco de la nave
$f_{sc} \sim 0.01\text{-}0.05$ Hz. Su sello es la **inversión de polarización** —
zurda (left-handed) en el marco de la nave pero intrínsecamente diestra
(right-handed) en el marco del plasma, un efecto Doppler porque el viento
super-alfvénico arrastra hacia atrás las ondas que se propagan hacia el Sol. La
genera la inestabilidad de haz resonante ion–ion de mano derecha
($\omega - \mathbf{k}\cdot\mathbf{V}_b = +\Omega_{cp}$). Otras poblaciones son las
ondas de Alfvén/ciclotrón iónico de "3 segundos" ($f_{sc} \sim 0.3$ Hz) y los
*whistlers* de "1 Hz" ($f_{sc} \sim 0.5\text{-}1$ Hz, diestros, generados por el
choque). Las ondas ULF pueden empinarse hasta formar *shocklets* y SLAMS
($\delta B/B_0 > 2$).

Corriente abajo, la vaina comprimida desarrolla $T_\perp > T_\parallel$, lo que
excita **ondas de Alfvén/ciclotrón iónico (AIC)** ($\beta_p < 1$; transversales,
zurdas, $\mathbf{k}\parallel\mathbf{B}$) y **modos espejo (mirror modes)**
($\beta_p > 1$; no propagantes, de polarización lineal, compresivos, con
$\delta B$ anticorrelacionado con $\delta n$, $\mathbf{k}$ casi-perpendicular).
Para el análisis espectral el contraste es nítido: el espectro **corriente
arriba** muestra **picos discretos** con fluctuaciones transversales
circulares/elípticas, mientras que la vaina **corriente abajo** es **de banda
ancha, más empinada y más compresiva** (más elíptica-a-lineal). Cualquier
inferencia del *sentido* de polarización debe corregirse por Doppler
($\omega_{sc} = \omega_{pf} + \mathbf{k}\cdot\mathbf{V}_{sw}$).

## (c) FFT y densidad espectral de potencia

Un magnetómetro entrega $B(t_n)$, $t_n = n\,\delta t$, $n = 0\dots N-1$, con
frecuencia de muestreo $f_s = 1/\delta t$ durante $T = N\,\delta t$. La
transformada de Fourier discreta (DFT)
$\tilde{B}_k = \sum_n B_n\, e^{-i 2\pi k n/N}$ lleva $N$ muestras a $N$
coeficientes complejos en $f_k = k/T$. La
**frecuencia de Nyquist** es $f_{\mathrm{Nyq}} = f_s/2$ (el contenido por encima
se *aliasa*); el espaciamiento es $\Delta f = 1/T = f_s/N$ — la resolución en
frecuencia solo mejora con un **registro más largo**, no con un muestreo más
fino. Una DFT directa cuesta $O(N^2)$; la **FFT de Cooley–Tukey** divide
recursivamente subsucesiones par/impar ($\tilde{B}_k = E_k + W_N^k O_k$) hasta
$O(N\log_2 N)$ — unas 102× más rápida con $N = 1024$, el motor de cada espectro
aquí.

El **teorema de Parseval** $\sum_n |B_n|^2 = \tfrac{1}{N}\sum_k |\tilde{B}_k|^2$
fija la normalización de la DEP para que el espectro porte la misma potencia que
la serie temporal. La DEP de un solo lado en nT²/Hz es

```math
S(f) = \frac{2}{f_s N}\,|\tilde{B}_k|^2 = \frac{2\,\delta t}{N}\,|\tilde{B}_k|^2 ,
```

donde el factor 2 pliega la potencia de frecuencias negativas sobre las positivas
(las celdas DC y de Nyquist conservan peso 1). Para una ventana $w_n$ se divide
entre $\sum_n w_n^2$ en vez de $N$. La varianza se recupera entonces como
$\sigma_B^2 = \sum_k S(f_k)\,\Delta f$, y la **DEP de la traza**
$S_{\mathrm{tr}} = S_x + S_y + S_z$ es la cantidad estándar en turbulencia. Un
registro finito convoluciona el espectro verdadero con un sinc cuyo primer lóbulo
lateral está solo $-13$ dB por debajo, produciendo **fuga espectral** y sesgando
los espectros de ley de potencia empinados; las **ventanas** suaves lo suprimen —
Hann $w(n) = 0.5 - 0.5\cos(2\pi n/N)$ ($-31$ dB), Hamming ($-43$ dB). El
**método de Welch** divide el registro en $K$ segmentos enventanados solapados
(típicamente 50 %) y promedia sus periodogramas, reduciendo la varianza del
estimador en $\sim 1/K$ a costa de una resolución más burda
$\Delta f = 1/(M\,\delta t)$ — el compromiso fundamental
**resolución ↔ estabilidad**.

## (d) Transformada wavelet continua y cono de influencia

Siguiendo a **Torrence & Compo (1998)**, la wavelet madre **Morlet** compleja es
una onda plana modulada por una gaussiana

```math
\psi_0(\eta) = \pi^{-1/4}\,e^{i\omega_0\eta}\,e^{-\eta^2/2},
```

usualmente con $\omega_0 = 6$. Al ser compleja entrega **amplitud y fase**. La
CWT es la convolución
$W_n(s) = \sum_{n'} x_{n'}\,\psi^{\ast}\!\big[(n'-n)\,\delta t/s\big]$, evaluada
eficientemente en el dominio de Fourier con cada escala normalizada a energía
unitaria. Las escalas siguen una malla diádica sub-octava
$s_j = s_0\,2^{j\,\delta j}$, escala mínima $s_0 = 2\,\delta t$ y
$\delta j = 1/8\dots 1/12$ para un espectro suave. La escala wavelet **no** es el
periodo de Fourier; para $\omega_0 = 6$ el periodo de Fourier equivalente es
$\lambda = 1.033\,s$ ($\approx$ escala), de modo que $f \approx \omega_0/(2\pi s)$
con un error de pocos por ciento.

El **cono de influencia (COI)** es la región donde los efectos de borde de la
serie finita (rellenada con ceros) corrompen la transformada; para Morlet está
acotado por el tiempo de caída $\tau_s = \sqrt{2}\,s$, se ensancha con la escala y
se **enmascara** porque la potencia de baja frecuencia allí no es confiable. La
DEP tiempo–frecuencia normalizada físicamente es

```math
S(s,t) = 2\,\delta t\,\mathrm{Tr}\{S_{ij}(f,t)\}
       = 2\,\delta t\,\big(|W_x|^2 + |W_y|^2 + |W_z|^2\big)\quad [\mathrm{nT^2/Hz}],
```

(más precisamente el prefactor
$\delta t\,F_F/(C_\delta\ln 2)\approx 1.92\,\delta t$, con $C_\delta = 0.776$). El
**promedio temporal** da el *espectro wavelet global* — un estimador insesgado y
**más suave** de la misma DEP que la FFT. La CWT respeta el límite de Heisenberg
$\delta f\,\delta t \ge 1/4\pi$ con un **embaldosado adaptativo ($Q$ constante)**:
buena resolución temporal a $f$ alta, buena resolución en frecuencia a $f$ baja —
superior a la FFT enventanada para señales no estacionarias como el cruce de un
choque.

## (e) Análisis de varianza mínima y polarización

El **MVA** (Sonnerup & Scheible 1998) estima la dirección normal / de propagación
de una estructura 1-D o de una onda plana como la dirección de **menor variación
del campo**. A partir de la matriz de varianza (covarianza)

```math
M_{ij} = \langle B_i B_j\rangle - \langle B_i\rangle\langle B_j\rangle ,
```

el problema de valores propios $M\,\mathbf{x} = \lambda\,\mathbf{x}$ da
$\lambda_1 \ge \lambda_2 \ge \lambda_3$ (máximo, intermedio, mínimo) con vectores
propios ortonormales: $\mathbf{x}_3$ (varianza mínima) $\approx$ el vector de onda
$\mathbf{k}$ / normal $\mathbf{n}$ de la estructura, mientras
$\mathbf{x}_1,\mathbf{x}_2$ generan el plano de polarización.

Diagnósticos a partir de las **razones de valores propios**:

| razón | significado |
|---|---|
| $\lambda_2/\lambda_3 \gg 1$ | la dirección de varianza mínima (plano de polarización / normal) está **bien definida** |
| $\lambda_1/\lambda_2 \approx 1$ | polarización **circular** |
| $\lambda_1/\lambda_2 \gg 1$ | polarización **lineal** |
| intermedio | **elíptica**, con elipticidad $\varepsilon \approx \sqrt{\lambda_2/\lambda_1}$ |

El **hodograma** ($\mathbf{B}$ proyectada en el plano
$\mathbf{x}_1$–$\mathbf{x}_2$) muestra la elipse directamente, y su sentido de
giro en torno a $\mathbf{k}$ da la quiralidad. Khrabrov & Sonnerup (1998) dan el
error angular

```math
|\Delta\varphi_{ij}| = \sqrt{\frac{\lambda_3\,(\lambda_i + \lambda_j - \lambda_3)}{(M-1)\,(\lambda_i - \lambda_j)^2}},
```

que diverge cuando $\lambda_i \to \lambda_3$ y decrece como $1/\sqrt{M}$; un
umbral práctico de confiabilidad es $\lambda_2/\lambda_3 \gtrsim 5\text{-}10$.

El **complemento resuelto en frecuencia** es la matriz espectral hermitiana
$S_{ij}(f) = \langle B_i(f)\,B_j^{\ast}(f)\rangle$ (aquí construida desde la CWT). De
ella se obtiene, por frecuencia, el **grado de polarización**

```math
\mathrm{DOP} = \sqrt{\frac{3\,\mathrm{Tr}S^2 - (\mathrm{Tr}S)^2}{2\,(\mathrm{Tr}S)^2}},
```

la **elipticidad** $2\,|\mathrm{Im}\,S|_{\mathrm{off}}/\mathrm{Tr}S$
(0 = lineal, 1 = circular; ambos invariantes rotacionales y por tanto
independientes del marco) y — vía la parte imaginaria (antisimétrica) de $S$,
Means (1972) — la **dirección normal de la onda**.

## (f) El conjunto de datos Wind / MFI

La nave **Wind** de la NASA (lanzada el 1994-11-01) es un monitor del viento solar
que mantiene una órbita de halo alrededor del punto **L1** Sol–Tierra
($\sim 1.5\times 10^6$ km hacia el Sol) desde 2004, muestreando el viento solar
prístino corriente arriba del choque de proa terrestre. El instrumento **MFI**
(Lepping et al. 1995) es un magnetómetro fluxgate triaxial doble (sensor externo
en una pértiga de 12 m) con 8 rangos dinámicos ($\pm 0.001$ a $\pm 65\,536$ nT),
digitización de 12 bits y ruido $< 0.006$ nT rms. El producto usado aquí,
**`WI_H2_MFI`** ("alta resolución"), tiene cadencia $\delta t \approx 0.092$ s
($\sim 11$ vec/s, **~39 000 puntos/hora**, Nyquist $\approx 5.4$ Hz) —
exactamente la cifra citada en las notas del curso. Las variables CDF son
`BGSE/BGSM` (campo vectorial de 3 componentes, nT) y `BF1` (escalar $|B|$). Las
coordenadas **GSE** ubican $+X$ hacia el Sol, $+Z$ al norte eclíptico,
$+Y = Z\times X$.

Los parámetros del choque provienen de la **Base de Datos de Choques
Interplanetarios del CfA** mediante el método **RH08** (Koval & Szabo 2008), un
ajuste por mínimos cuadrados no lineales de las ocho relaciones de
Rankine–Hugoniot para $\hat{\mathbf{n}}$ y $V_{sh}$. Para el **CfA #00802**
(2022-08-19 16:51 UT): tipo FF; $\hat{\mathbf{n}} = (-0.994,\,0.114,\,0.006)$ GSE;
$\theta_{Bn} = 31.3 \pm 11.6^\circ$; $V_{sh} = 683.4 \pm 2.7$ km/s;
$r = 2.39 \pm 0.24$; corriente arriba $\mathbf{B} = (4.7,\,1.6,\,-1.1)$ nT,
$n = 3.1\ \mathrm{cm^{-3}}$, $v_A \approx 67$ km/s, $c_s \approx 64$ km/s,
$M_f = 1.81 \pm 0.08$; corriente abajo $\mathbf{B} = (5.4,\,2.8,\,-7.8)$ nT
($B_z$ fuertemente negativo, es decir geoefectivo), $n = 7.3\ \mathrm{cm^{-3}}$.
Los datos los distribuye NASA SPDF / CDAWeb como archivos CDF autodescriptivos
según las convenciones ISTP.

---

## Referencias

- Treumann (2009), *Fundamentals of collisionless shocks for astrophysical application I*, A&ARv 17, 409. <https://link.springer.com/article/10.1007/s00159-009-0024-2>
- Burgess (2008), *Non-relativistic Collisionless Shock Physics IV: Quasi-Parallel Supercritical Shocks*. <https://arxiv.org/abs/0805.2579>
- Kilpua et al. (2015), *Properties and drivers of fast interplanetary shocks near Earth (1995–2013)*, JGR 120. <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JA021138>
- Eastwood et al. (2005), *The Foreshock*, Space Sci. Rev. <https://space.physics.uiowa.edu/~dag/publications/2005_TheForeshock_SSR.pdf>
- Hoppe & Russell (1983), *Plasma rest frame frequencies and polarizations of low-frequency upstream waves*, JGR. <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JA088iA03p02021>
- Paschmann & Daly (eds.) (1998), *Analysis Methods for Multi-Spacecraft Data*, ISSI SR-001. <https://www.issibern.ch/PDF-Files/analysis_methods_1_1a.pdf>
- Torrence & Compo (1998), *A Practical Guide to Wavelet Analysis*, BAMS 79(1), 61–78. <https://psl.noaa.gov/people/gilbert.p.compo/Torrence_compo1998.pdf>
- Khrabrov & Sonnerup (1998), *Error estimates for minimum variance analysis*, JGR 103(A4). <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/97JA03731>
- Santolík, Parrot & Lefeuvre (2003), *SVD methods for wave propagation analysis*, Radio Sci. 38(1). <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2000RS002523>
- Lepping et al. (1995), *The Wind Magnetic Field Investigation*, Space Sci. Rev. 71, 207. <https://wind.nasa.gov/docs/MFI_Lepping_SSR1995.pdf>
- Koval & Szabo (2008), *Modified Rankine–Hugoniot shock fitting technique (RH08)*, JGR. <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JA013337>
- Base de Datos de Choques IP del CfA — Wind #00802. <https://lweb.cfa.harvard.edu/shocks/wi_data/00802/wi_00802.html>
- Algoritmo FFT de Cooley–Tukey. <https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm>
- Método de Welch. <https://en.wikipedia.org/wiki/Welch%27s_method>
