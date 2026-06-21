<div align="center">

# 🌌 Spectral Analysis of an Interplanetary Shock
### Análisis Espectral de un Choque Interplanetario

**FFT · Continuous Wavelet Transform · Minimum Variance Analysis**
of high-resolution *Wind*/MFI magnetic-field data across a fast-forward shock.

<br>

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![NumPy](https://img.shields.io/badge/NumPy-013243?style=for-the-badge&logo=numpy&logoColor=white)](https://numpy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)](https://scipy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?style=for-the-badge)](https://matplotlib.org/)
[![License](https://img.shields.io/badge/License-MIT-3DA639?style=for-the-badge)](https://opensource.org/license/mit)

[![Data](https://img.shields.io/badge/Data-Wind%2FMFI%20·%20CDAWeb-FF6F00?style=flat-square)](https://cdaweb.gsfc.nasa.gov/)
[![Event](https://img.shields.io/badge/shock-2022--08--19%2016%3A51%20UT-8957e5?style=flat-square)](https://lweb.cfa.harvard.edu/shocks/wi_data/00802/wi_00802.html)
[![Methods](https://img.shields.io/badge/methods-FFT%20·%20CWT%20·%20MVA-1f6feb?style=flat-square)](docs/theory_en.md)
[![Tests](https://img.shields.io/badge/tests-6%20passing-2ea44f?style=flat-square)](tests/test_analysis.py)
[![Reproducible](https://img.shields.io/badge/reproducible-offline%20cache-success?style=flat-square)](scripts/run_analysis.py)
[![Course](https://img.shields.io/badge/UNAM-Instituto%20de%20Geofísica-9d2449?style=flat-square)](https://www.geofisica.unam.mx/)

<br>

![Overview of the shock](figures/fig01_overview.png)

<sub>NASA *Wind*/MFI magnetic field across a fast-forward interplanetary shock —
quiet, wave-filled **upstream** (blue) vs. compressed, turbulent **downstream** (red).</sub>

**🌐 [Español](README.md)  ·  English (this file)**

</div>

---

## Overview

This repository implements the homework of the UNAM course *Space Science and
Data Analysis — 3. Spectral Analysis* (Dr. Byeongseon Park, Instituto de
Geofísica, UNAM). For one real interplanetary shock it:

1. picks a **shock event** and defines **upstream** `[-11, -1] min` and
   **downstream** `[+1, +11] min` intervals around it;
2. computes the **power spectral density (PSD)** of the magnetic field with both
   the **FFT** (periodogram + Welch) and a **Morlet continuous wavelet transform
   (CWT)**;
3. produces the **2-D time–frequency** wave-activity maps (CWT scalograms);
4. analyses the **wave polarization** of each interval with **minimum variance
   analysis (MVA)** and a frequency-resolved spectral matrix.

> 📄 Full scientific background: [`docs/theory_en.md`](docs/theory_en.md) ·
> assignment: [`docs/space_science_data_analysis_3.pdf`](docs/space_science_data_analysis_3.pdf)

## The event

A **fast-forward (FF)** interplanetary shock observed by **Wind** at the L1
Lagrange point, taken from the [CfA Interplanetary Shock Database
(#00802)](https://lweb.cfa.harvard.edu/shocks/wi_data/00802/wi_00802.html).

| Quantity | Value | Meaning |
|---|---|---|
| Time `t₀` | **2022-08-19 16:51:00 UT** | shock arrival at Wind |
| Type | **FF** | fast-forward (B and n increase) |
| `θ_Bn` | **31.3° ± 11.6°** | **quasi-parallel** geometry |
| `V_shock` | 683 km/s | shock speed |
| `M_fast` | 1.81 | fast magnetosonic Mach number |
| `β_up` | 0.55 | upstream plasma beta |
| `r = B_d/B_u` | ≈ 2.4 (CfA) / ≈ 2.0 (obs.) | magnetic compression |
| Dataset | `WI_H2_MFI` (`BGSE`) | high-res field, **fs ≈ 10.87 Hz** (≈ 39 000 pts/h) |

The data are fetched on-the-fly from **NASA/GSFC CDAWeb** (a small subset CDF),
read with `cdflib`, and cached as a compressed `.npz` so every figure
**reproduces offline**.

## Methods in one paragraph

The **FFT PSD** uses the one-sided convention
$S(f) = 2\,\delta t/N\,|\tilde{B}(f)|^2$ in nT²/Hz, normalized so that
$\int S\,df = \mathrm{var}(B)$ (verified to **1.000**). The **CWT** follows
*Torrence & Compo (1998)*: complex Morlet wavelet
$\psi_0(\eta)=\pi^{-1/4}e^{i\omega_0\eta}e^{-\eta^2/2}$ ($\omega_0=6$), dyadic
scales $s_j=s_0\,2^{j\,\delta j}$, cone of influence $\tau_s=\sqrt{2}\,s$, and the
time–frequency PSD $S(f,t)=2\,\delta t\,\mathrm{Tr}\{S_{ij}\}$ with
$S_{ij}=W_i W_j^{\ast}$. **MVA** diagonalizes the variance matrix
$M_{ij}=\langle B_iB_j\rangle-\langle B_i\rangle\langle B_j\rangle$; the eigenvalue
ratios $\lambda_1/\lambda_2$ (linear↔circular) and $\lambda_2/\lambda_3$ (quality
of the normal) classify the wave polarization.

## Results

### 1 · Power spectral density — FFT vs CWT

![PSD FFT vs CWT](figures/fig02_psd_fft_cwt.png)

The Welch-averaged FFT (color) and the global wavelet spectrum (green) **agree**,
but the **CWT is markedly smoother** — the key point of the lecture. The
**downstream** spectrum carries ~4× more power and shows a developed turbulent
cascade close to the **Kolmogorov $f^{-5/3}$** reference; the **upstream**
spectrum is flatter with a ULF foreshock shoulder near 0.02–0.06 Hz.

### 2 · Time–frequency wave activity (CWT scalograms)

![CWT scalograms](figures/fig03_cwt_scalogram.png)

Power jumps sharply at the shock (`t = 0`). Upstream the power is concentrated at
low frequencies (the ~30-s foreshock ULF waves); downstream it becomes
**broadband**, extending to higher frequencies (the turbulent sheath). The
hatched **cone of influence** marks edge-affected (unreliable) regions.

### 3 · Minimum variance analysis (hodograms)

![MVA hodograms](figures/fig04_mva_hodograms.png)

Band-passed to the wave band (0.02–0.25 Hz). **Upstream** is nearly **circular**
(`λ₁/λ₂ = 1.30`, well-defined plane `λ₂/λ₃ = 5.1`, propagating ~along **B**,
θ_kB ≈ 11°) — the signature of quasi-parallel **foreshock ULF waves**.
**Downstream** is clearly **elliptical** (`λ₁/λ₂ = 2.13`) and larger amplitude —
compressive sheath turbulence. This reproduces the lecture's slide-11 example
(upstream circular, downstream elliptical).

### 4 · Frequency-resolved polarization *(bonus, slide 9)*

![Polarization](figures/fig05_polarization.png)

Degree of polarization and ellipticity from the CWT spectral matrix `S_ij(f,t)`
(time-averaged over ~5 wave periods per scale), a frequency-resolved complement
to MVA.

### Key numbers

| Interval | window | ⟨\|B\|⟩ | var | λ₁/λ₂ | λ₂/λ₃ | ellipticity | polarization | θ_kB |
|---|---|---|---|---|---|---|---|---|
| **Upstream** | [−11,−1] min | 5.52 nT | 5.2 nT² | **1.30** | 5.08 | 0.88 | **circular** | 11° |
| **Downstream** | [+1,+11] min | 10.99 nT | 20.1 nT² | **2.13** | 2.64 | 0.69 | **elliptical** | 4° |

<sub>Parseval check: FFT ∫S df / var = 1.000 · CWT reconstruction / var = 0.90
(expected Torrence–Compo bias).</sub>

## Quick start

```bash
# 1. set up the environment
make setup                     # or: python -m venv .venv && .venv/bin/pip install -r requirements.txt

# 2. run the whole pipeline (downloads + caches data, writes figures/*.png)
make run                       # or: .venv/bin/python scripts/run_analysis.py

# offline / no network? use the physical synthetic fallback:
make synthetic

# run the tests
make test
```

Analyse a different shock by editing the `EVENT` object in
[`src/spectral_shock/config.py`](src/spectral_shock/config.py).

## Repository structure

```
.
├── README.md · README.en.md      # Spanish (default) + this English version
├── requirements.txt · Makefile · LICENSE
├── scripts/run_analysis.py       # end-to-end pipeline → figures/ + results.json
├── src/spectral_shock/
│   ├── config.py                 # event, intervals, parameters
│   ├── data.py                   # CDAWeb download + cdflib + cache + fallback
│   ├── spectral.py               # FFT PSD, Morlet CWT, COI, spectral matrix
│   ├── mva.py                    # minimum variance analysis
│   └── plotting.py               # all figures
├── tests/test_analysis.py        # offline unit tests (Parseval, MVA, …)
├── data/                         # cached .npz (committed) + .cdf (ignored)
├── figures/                      # generated PNGs + results.md
└── docs/                         # theory_en.md, theory_es.md, assignment PDF
```

---

<div align="center">

### 📚 References · 🛰️ Data

</div>

**Methods:** Torrence & Compo (1998), *A Practical Guide to Wavelet Analysis*,
BAMS; Sonnerup & Scheible (1998) in Paschmann & Daly (eds.), *Analysis Methods
for Multi-Spacecraft Data*, ISSI SR-001; Khrabrov & Sonnerup (1998). Full lists
in [`docs/theory_en.md`](docs/theory_en.md) /
[`docs/theory_es.md`](docs/theory_es.md).

**Data:** NASA *Wind*/MFI magnetic field (Lepping et al. 1995), dataset
`WI_H2_MFI`, courtesy of the Wind MFI team and **NASA/GSFC SPDF
[CDAWeb](https://cdaweb.gsfc.nasa.gov/)**. Shock parameters from the
**[CfA Interplanetary Shock Database](https://lweb.cfa.harvard.edu/shocks/)**
(RH08, Koval & Szabo 2008). Please cite these providers when using the data.

**License:** Code under the [MIT License](LICENSE). NASA/CfA data are public and
not covered by this license.

<div align="center">
<sub>Built for the UNAM <i>Space Science and Data Analysis</i> course · Author: Alonso Cervantes Flores</sub>
</div>
