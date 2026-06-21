# Scientific background — spectral analysis of a collisionless shock

> Companion notes for the repository **Spectral analysis of an interplanetary
> shock**. They expand the condensed summary in the main `README.md`.
> A Spanish version is available in [`theory_es.md`](theory_es.md).
>
> *(Equations are written in GitHub-flavoured LaTeX and render on github.com.)*

---

## (a) Collisionless and interplanetary shocks

A shock is an irreversible transition that decelerates a super-magnetosonic flow
to sub-magnetosonic speed while **compressing and heating** the plasma. In an
ordinary gas the dissipation is supplied by Coulomb collisions and the front is a
few mean free paths thick. In the solar wind the Coulomb mean free path is of
order **1 AU** ($\sim 10^8$ km), yet observed shock fronts are only
$\sim 10^2\text{-}10^3$ km thick. The dissipation must therefore be
**collisionless** — provided by collective electromagnetic processes:
wave–particle interactions, micro-instabilities, and ion reflection (an
"anomalous" resistivity/viscosity) that thermalise the flow on kinetic scales.
The relevant lengths are the ion inertial length $l_i = c/\omega_{pi}$
($\approx 228/\sqrt{n\,[\mathrm{cm^{-3}}]}$ km, $\approx 130$ km at
$n = 3.1\ \mathrm{cm^{-3}}$) and the ion gyroradius $r_{ci} = v_{th}/\Omega_{ci}$.

Ideal MHD admits three wave modes: the compressive **fast** and **slow**
magnetosonic modes and the non-compressive intermediate **Alfvén** mode
($v_A = B/\sqrt{\mu_0\rho}$). The magnetosonic phase speeds are

```math
v_{f,s}^2 = \tfrac{1}{2}\Big[(v_A^2 + c_s^2) \pm \sqrt{(v_A^2 + c_s^2)^2 - 4\,v_A^2 c_s^2 \cos^2\theta}\,\Big],
\qquad c_s = \sqrt{\gamma p/\rho}.
```

Interplanetary (IP) shocks are overwhelmingly **fast magnetosonic** shocks. A
**fast-forward (FF)** shock propagates anti-sunward faster than the upstream wind
(the most common type at 1 AU, usually driven by ICMEs or stream interaction
regions); **B and density increase across the front** — exactly our 2022-08-19
event. The **Rankine–Hugoniot** relations conserve mass, momentum, energy,
$\nabla\cdot\mathbf{B}=0$ and the frozen-in condition across the front. In the
gas-dynamic limit the compression is
$r = \rho_d/\rho_u = (\gamma+1)M^2/[(\gamma-1)M^2 + 2]$, with strong-shock maximum
$r = (\gamma+1)/(\gamma-1) = 4$ for $\gamma = 5/3$; the magnetic compression
$B_d/B_u$ is likewise $\le 4$. Real 1 AU shocks have $r \approx 1.5\text{-}3$ —
our event has a **moderate** $r \approx 2.4$ (CfA fit), consistent with the
observed $|B|$ jump from $\approx 5.5$ to $\approx 11$ nT.

The **shock-normal angle** $\theta_{Bn}$ (between the upstream $\mathbf{B}$ and
the shock normal $\mathbf{n}$) splits shocks into **quasi-perpendicular**
($\theta_{Bn} > 45^\circ$; thin, laminar, foot–ramp–overshoot) and
**quasi-parallel** ($\theta_{Bn} < 45^\circ$; extended, turbulent, with reflected
ions streaming far upstream to form a *foreshock* and SLAMS). At
$\theta_{Bn} \approx 31^\circ$ our shock is **quasi-parallel**. A shock requires
fast Mach number $M_f = v_{un}/v_f > 1$; the **first critical Mach number**
($1 \le M_c \le 2.76$) separates subcritical (resistive) from supercritical
(ion-reflection-dominated) shocks. With $M_f \approx 1.81$ this is a
weak-to-moderate, near-critical quasi-parallel FF shock.

## (b) Upstream vs downstream wave activity

Geometry sets the wave physics. Quasi-parallel shocks let ions escape sunward
along the field, producing an **extended ion foreshock** (field-aligned beams →
intermediate → gyrating → diffuse ions) that is filled with waves. The dominant
mode is the **"30-second" ULF wave**: quasi-monochromatic, large-amplitude
($\delta B/B \sim 0.1\text{-}1$), fast-magnetosonic, at spacecraft-frame
frequency $f_{sc} \sim 0.01\text{-}0.05$ Hz. Its hallmark is **polarization
reversal** — left-handed in the spacecraft frame but intrinsically right-handed
in the plasma frame, a Doppler effect because the super-Alfvénic wind sweeps
sunward-propagating waves back past the spacecraft. It is driven by the ion–ion
right-hand resonant beam instability
($\omega - \mathbf{k}\cdot\mathbf{V}_b = +\Omega_{cp}$). Additional populations
are "3-second" Alfvén/ion-cyclotron waves
($f_{sc} \sim 0.3$ Hz) and "1-Hz" whistlers ($f_{sc} \sim 0.5\text{-}1$ Hz,
right-handed, shock-generated). ULF waves can steepen into shocklets and SLAMS
($\delta B/B_0 > 2$).

Downstream, the compressed sheath develops $T_\perp > T_\parallel$, driving
**Alfvén/ion-cyclotron (AIC) waves** ($\beta_p < 1$; transverse, left-handed,
$\mathbf{k}\parallel\mathbf{B}$) and **mirror modes** ($\beta_p > 1$;
non-propagating, linearly polarized, compressive, $\delta B$ anti-correlated with
$\delta n$, $\mathbf{k}$ quasi-perpendicular). For spectral analysis the contrast
is sharp: the **upstream** spectrum shows **discrete peaks** with transverse,
circular/elliptical fluctuations, whereas the **downstream** sheath is
**broadband, steeper and more compressive** (more elliptical-to-linear). Any
inference of the polarization *sense* must be Doppler-corrected
($\omega_{sc} = \omega_{pf} + \mathbf{k}\cdot\mathbf{V}_{sw}$).

## (c) FFT and Power Spectral Density

A magnetometer returns $B(t_n)$, $t_n = n\,\delta t$, $n = 0\dots N-1$, with
sampling rate $f_s = 1/\delta t$ over $T = N\,\delta t$. The discrete Fourier
transform (DFT) $\tilde{B}_k = \sum_n B_n\, e^{-i 2\pi k n/N}$ maps $N$ samples
to $N$ complex coefficients at $f_k = k/T$. The **Nyquist frequency** is
$f_{\mathrm{Nyq}} = f_s/2$ (content above it aliases); the bin spacing is
$\Delta f = 1/T = f_s/N$ — frequency resolution improves only with a **longer
record**, not finer sampling. A direct DFT costs $O(N^2)$; the **Cooley–Tukey
FFT** recursively splits even/odd subsequences ($\tilde{B}_k = E_k + W_N^k O_k$)
for $O(N\log_2 N)$ — about 102× faster at $N = 1024$, the engine behind every
spectrum here.

**Parseval's theorem** $\sum_n |B_n|^2 = \tfrac{1}{N}\sum_k |\tilde{B}_k|^2$ fixes
the PSD normalization so the spectrum carries the same power as the time series.
The one-sided PSD in nT²/Hz is

```math
S(f) = \frac{2}{f_s N}\,|\tilde{B}_k|^2 = \frac{2\,\delta t}{N}\,|\tilde{B}_k|^2 ,
```

where the factor 2 folds negative-frequency power onto positive frequencies (the
DC and Nyquist bins keep weight 1). For a window $w_n$, divide by $\sum_n w_n^2$
instead of $N$. The variance is then recovered as
$\sigma_B^2 = \sum_k S(f_k)\,\Delta f$, and the **trace PSD**
$S_{\mathrm{tr}} = S_x + S_y + S_z$ is the standard turbulence quantity. A finite
record convolves the true spectrum with a sinc whose first sidelobe is only
$-13$ dB down, **leaking** power and biasing steep power-law spectra; smooth
**windows** suppress this — Hann $w(n) = 0.5 - 0.5\cos(2\pi n/N)$ ($-31$ dB),
Hamming ($-43$ dB). **Welch's method** splits the record into $K$ overlapping
(typically 50 %) windowed segments and averages their periodograms, reducing the
estimator variance by $\sim 1/K$ at the cost of coarser resolution
$\Delta f = 1/(M\,\delta t)$ — the fundamental **resolution ↔ stability**
trade-off.

## (d) Continuous Wavelet Transform and the cone of influence

Following **Torrence & Compo (1998)**, the complex **Morlet** mother wavelet is a
Gaussian-modulated plane wave

```math
\psi_0(\eta) = \pi^{-1/4}\,e^{i\omega_0\eta}\,e^{-\eta^2/2},
```

commonly with $\omega_0 = 6$. Being complex it yields **amplitude and phase**.
The CWT is the convolution
$W_n(s) = \sum_{n'} x_{n'}\,\psi^{\ast}\!\big[(n'-n)\,\delta t/s\big]$, evaluated
efficiently in the Fourier domain with each scale normalized to unit energy.
Scales follow a sub-octave dyadic grid $s_j = s_0\,2^{j\,\delta j}$, smallest
scale $s_0 = 2\,\delta t$ and $\delta j = 1/8\dots 1/12$ for a smooth spectrum.
The wavelet scale is **not** the Fourier period; for $\omega_0 = 6$ the equivalent
Fourier period is $\lambda = 1.033\,s$ ($\approx$ scale), so
$f \approx \omega_0/(2\pi s)$ to within a few percent.

The **cone of influence (COI)** is the region where edge effects from the finite,
zero-padded series corrupt the transform; for Morlet it is bounded by the
e-folding time $\tau_s = \sqrt{2}\,s$, widens with scale, and is **masked**
because the low-frequency power there is unreliable. The physically normalized
time–frequency PSD is

```math
S(s,t) = 2\,\delta t\,\operatorname{Tr}\{S_{ij}(f,t)\}
       = 2\,\delta t\,\big(|W_x|^2 + |W_y|^2 + |W_z|^2\big)\quad [\mathrm{nT^2/Hz}],
```

(more precisely the prefactor $\delta t\,F_F/(C_\delta\ln 2)\approx 1.92\,\delta t$,
with $C_\delta = 0.776$). **Time-averaging** gives the *global wavelet spectrum* —
an unbiased, **smoother** estimator of the same PSD as the FFT. The CWT respects
the Heisenberg limit $\delta f\,\delta t \ge 1/4\pi$ with **adaptive
(constant-$Q$) tiling**: good time resolution at high $f$, good frequency
resolution at low $f$ — superior to the windowed FT for non-stationary signals
such as a shock crossing.

## (e) Minimum Variance Analysis and polarization

**MVA** (Sonnerup & Scheible 1998) estimates the normal / propagation direction
of a 1-D structure or planar wave as the direction of **least field variation**.
From the variance (covariance) matrix

```math
M_{ij} = \langle B_i B_j\rangle - \langle B_i\rangle\langle B_j\rangle ,
```

the eigenproblem $M\,\mathbf{x} = \lambda\,\mathbf{x}$ gives
$\lambda_1 \ge \lambda_2 \ge \lambda_3$ (maximum, intermediate, minimum) with
orthonormal eigenvectors: $\mathbf{x}_3$ (minimum-variance) $\approx$ the wave
vector $\mathbf{k}$ / structure normal $\mathbf{n}$, while
$\mathbf{x}_1,\mathbf{x}_2$ span the polarization plane.

Diagnostics from the **eigenvalue ratios**:

| ratio | meaning |
|---|---|
| $\lambda_2/\lambda_3 \gg 1$ | minimum-variance direction (polarization plane / normal) **well defined** |
| $\lambda_1/\lambda_2 \approx 1$ | **circular** polarization |
| $\lambda_1/\lambda_2 \gg 1$ | **linear** polarization |
| in between | **elliptical**, with ellipticity $\varepsilon \approx \sqrt{\lambda_2/\lambda_1}$ |

The **hodogram** ($\mathbf{B}$ projected on the $\mathbf{x}_1$–$\mathbf{x}_2$
plane) shows the ellipse directly, and its rotation sense about $\mathbf{k}$ gives
the handedness. Khrabrov & Sonnerup (1998) give the angular error

```math
|\Delta\varphi_{ij}| = \sqrt{\frac{\lambda_3\,(\lambda_i + \lambda_j - \lambda_3)}{(M-1)\,(\lambda_i - \lambda_j)^2}},
```

which blows up as $\lambda_i \to \lambda_3$ and shrinks as $1/\sqrt{M}$; a
practical reliability threshold is $\lambda_2/\lambda_3 \gtrsim 5\text{-}10$.

The **frequency-resolved complement** is the Hermitian spectral matrix
$S_{ij}(f) = \langle B_i(f)\,B_j^{\ast}(f)\rangle$ (here built from the CWT). From it
one obtains, per frequency, the **degree of polarization**

```math
\mathrm{DOP} = \sqrt{\frac{3\,\operatorname{Tr}S^2 - (\operatorname{Tr}S)^2}{2\,(\operatorname{Tr}S)^2}},
```

the **ellipticity** $2\,|\mathrm{Im}\,S|_{\mathrm{off}}/\operatorname{Tr}S$
(0 = linear, 1 = circular; both rotational invariants, hence frame-independent),
and — via the imaginary (antisymmetric) part of $S$, Means (1972) — the
**wave-normal direction**.

## (f) The Wind / MFI dataset

NASA's **Wind** spacecraft (launched 1994-11-01) is a solar-wind monitor that has
held a halo orbit about the Sun–Earth **L1** point ($\sim 1.5\times 10^6$ km
sunward of Earth) since 2004, sampling the pristine solar wind upstream of Earth's
bow shock. The **MFI** instrument (Lepping et al. 1995) is a dual triaxial
fluxgate magnetometer (outer sensor on a 12 m boom) with 8 dynamic ranges
($\pm 0.001$ to $\pm 65\,536$ nT), 12-bit digitization and $< 0.006$ nT rms
noise. The product used here, **`WI_H2_MFI`** ("high-resolution"), has cadence
$\delta t \approx 0.092$ s ($\sim 11$ vec/s, **~39 000 points/hour**, Nyquist
$\approx 5.4$ Hz) — exactly the figure quoted in the lecture notes. CDF variables
are `BGSE/BGSM` (3-vector field, nT) and `BF1` (scalar $|B|$). **GSE** coordinates
place $+X$ toward the Sun, $+Z$ to ecliptic north, $+Y = Z\times X$.

The shock parameters come from the **CfA Interplanetary Shock Database** via the
**RH08** method (Koval & Szabo 2008), a nonlinear least-squares fit of the eight
Rankine–Hugoniot relations for $\hat{\mathbf{n}}$ and $V_{sh}$. For **CfA #00802**
(2022-08-19 16:51 UT): type FF; $\hat{\mathbf{n}} = (-0.994,\,0.114,\,0.006)$ GSE;
$\theta_{Bn} = 31.3 \pm 11.6^\circ$; $V_{sh} = 683.4 \pm 2.7$ km/s;
$r = 2.39 \pm 0.24$; upstream $\mathbf{B} = (4.7,\,1.6,\,-1.1)$ nT,
$n = 3.1\ \mathrm{cm^{-3}}$, $v_A \approx 67$ km/s, $c_s \approx 64$ km/s,
$M_f = 1.81 \pm 0.08$; downstream $\mathbf{B} = (5.4,\,2.8,\,-7.8)$ nT (strongly
negative $B_z$, i.e. geoeffective), $n = 7.3\ \mathrm{cm^{-3}}$. Data are
distributed by NASA SPDF / CDAWeb as self-describing CDF files following ISTP
conventions.

---

## References

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
- CfA Interplanetary Shock Database — Wind #00802. <https://lweb.cfa.harvard.edu/shocks/wi_data/00802/wi_00802.html>
- Cooley–Tukey FFT algorithm. <https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm>
- Welch's method. <https://en.wikipedia.org/wiki/Welch%27s_method>
