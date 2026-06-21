r"""
Spectral analysis
=================

Two estimators of the magnetic-field power spectral density (PSD), in
``nT^2 / Hz``, plus the time-frequency wavelet scalogram.

1. **FFT PSD** (Fourier).  One-sided periodogram and Welch-averaged spectrum,
   normalised so that :math:`\int S(f)\,df = \mathrm{var}(B)` (Parseval).
   This is the ``S(f) = 2T\,|\tilde B(n)|^2`` of the lecture notes.

2. **Morlet CWT** following Torrence & Compo (1998).  Complex Morlet mother
   wavelet

   .. math::  \psi_0(\eta) = \pi^{-1/4}\, e^{i\omega_0\eta}\, e^{-\eta^2/2},

   evaluated efficiently in the Fourier domain.  We build the scales
   :math:`s_j = s_0 2^{j\,\delta j}`, the cone of influence (COI), the
   3-component spectral matrix :math:`S_{ij}(f,t)=W_i W_j^{*}`, the scalogram
   :math:`S(f,t)=2\,\delta t\,\mathrm{Tr}\{S_{ij}\}` and the time-averaged
   (global) wavelet PSD, which is directly comparable with the FFT PSD.

The CWT PSD normalisation uses the bias-corrected constant of Torrence & Compo
(:math:`C_\delta = 0.776` for Morlet :math:`\omega_0 = 6`); the resulting
prefactor :math:`\delta t\,F_F/(C_\delta\ln 2)\simeq 1.92\,\delta t \approx
2\,\delta t`, matching the lecture's ``S = 2\delta t\,\mathrm{Tr}\{S_{ij}\}``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

# Morlet bias-correction constant for omega0 = 6 (Torrence & Compo 1998, Table 2)
C_DELTA = {6.0: 0.776}


# --------------------------------------------------------------------------- #
# FFT power spectral density
# --------------------------------------------------------------------------- #
@dataclass
class FFTSpectrum:
    f: np.ndarray                      # frequency [Hz]
    psd_comp: np.ndarray               # (3, Nf) per-component PSD [nT^2/Hz]
    psd_total: np.ndarray              # (Nf,) trace PSD = sum of components
    f_raw: np.ndarray                  # periodogram frequency
    psd_total_raw: np.ndarray          # raw (un-averaged) trace PSD
    method: str = "welch"


def fft_psd(
    B: np.ndarray,
    dt: float,
    *,
    seg_s: float = 100.0,
    overlap: float = 0.5,
    window: str = "hann",
    detrend: str = "linear",
) -> FFTSpectrum:
    """One-sided PSD of a (N,3) field via Welch + raw periodogram.

    ``scaling='density'`` returns nT^2/Hz and satisfies Parseval:
    ``trapz(psd, f) == var(B)`` (to within windowing bias).
    """
    B = np.asarray(B, dtype=float)
    B = _fill_gaps(B)
    n = B.shape[0]
    fs = 1.0 / dt
    nperseg = int(min(n, max(64, round(seg_s * fs))))
    noverlap = int(nperseg * overlap)

    psd_comp = []
    for c in range(3):
        f, p = signal.welch(B[:, c], fs=fs, window=window, nperseg=nperseg,
                            noverlap=noverlap, detrend=detrend, scaling="density")
        psd_comp.append(p)
    psd_comp = np.array(psd_comp)
    psd_total = psd_comp.sum(axis=0)

    # raw periodogram (boxcar, no averaging) of the trace power
    fr, pr = signal.periodogram(B, fs=fs, window="boxcar", detrend=detrend,
                                scaling="density", axis=0)
    psd_total_raw = pr.sum(axis=1)

    return FFTSpectrum(f=f, psd_comp=psd_comp, psd_total=psd_total,
                       f_raw=fr, psd_total_raw=psd_total_raw)


# --------------------------------------------------------------------------- #
# Morlet continuous wavelet transform (Torrence & Compo 1998)
# --------------------------------------------------------------------------- #
@dataclass
class CWT:
    W: np.ndarray            # (3, Ns, N) complex wavelet coefficients (per comp)
    scales: np.ndarray       # (Ns,) wavelet scales [s]
    freqs: np.ndarray        # (Ns,) Fourier-equivalent frequency [Hz]
    coi: np.ndarray          # (N,) cone-of-influence Fourier period [s]
    t: np.ndarray            # (N,) time [s]
    dt: float
    w0: float
    dj: float

    # --- scalogram & PSD --------------------------------------------------
    @property
    def trace_power(self) -> np.ndarray:
        """Tr{S_ij} = |Wx|^2+|Wy|^2+|Wz|^2  -> (Ns, N) [nT^2]."""
        return np.sum(np.abs(self.W) ** 2, axis=0)

    @property
    def scalogram(self) -> np.ndarray:
        """Time-frequency PSD  S(f,t) ~ 2 dt Tr{S_ij}  -> (Ns, N) [nT^2/Hz]."""
        return self._psd_prefactor() * self.trace_power

    def _psd_prefactor(self) -> float:
        ff = fourier_factor(self.w0)
        cdelta = C_DELTA.get(self.w0, 0.776)
        return self.dt * ff / (cdelta * np.log(2.0))

    def coi_mask(self) -> np.ndarray:
        """Boolean mask (Ns, N): True where INSIDE the COI (unreliable)."""
        period = 1.0 / self.freqs
        return period[:, None] > self.coi[None, :]

    def global_psd(self, mask_coi: bool = True) -> tuple[np.ndarray, np.ndarray]:
        """Time-averaged wavelet PSD over the record -> (freqs, psd) [nT^2/Hz]."""
        power = self.trace_power.copy()
        if mask_coi:
            m = self.coi_mask()
            power = np.where(m, np.nan, power)
        with np.errstate(invalid="ignore"):
            valid = np.sum(~np.isnan(power), axis=1)
            gp = np.where(valid > 0, np.nansum(np.nan_to_num(power), axis=1)
                          / np.maximum(valid, 1), np.nan)
        return self.freqs, self._psd_prefactor() * gp


def fourier_factor(w0: float) -> float:
    """Fourier period / scale ratio for the Morlet wavelet."""
    return (4.0 * np.pi) / (w0 + np.sqrt(2.0 + w0 ** 2))


def make_scales(n: int, dt: float, dj: float, s0_factor: float, w0: float):
    """Logarithmic scale set and the matching Fourier frequencies."""
    s0 = s0_factor * dt
    j_max = int(np.floor(np.log2(n * dt / s0) / dj))
    j = np.arange(j_max + 1)
    scales = s0 * 2.0 ** (j * dj)
    freqs = 1.0 / (fourier_factor(w0) * scales)
    return scales, freqs


def _morlet_ft(scale: float, omega: np.ndarray, dt: float, w0: float) -> np.ndarray:
    """Daughter wavelet in the Fourier domain (Torrence & Compo eq. 6)."""
    heav = omega > 0
    norm = np.sqrt(2.0 * np.pi * scale / dt) * (np.pi ** -0.25)
    expnt = -0.5 * (scale * omega - w0) ** 2
    return norm * np.exp(expnt) * heav


def cwt_component(x: np.ndarray, dt: float, scales: np.ndarray, w0: float) -> np.ndarray:
    """Morlet CWT of a single 1-D series -> (Ns, N) complex."""
    x = np.asarray(x, dtype=float)
    x = _fill_gaps_1d(x)
    x = x - np.mean(x)
    n = x.size
    xhat = np.fft.fft(x)
    omega = 2.0 * np.pi * np.fft.fftfreq(n, d=dt)
    W = np.empty((scales.size, n), dtype=complex)
    for i, s in enumerate(scales):
        daughter = _morlet_ft(s, omega, dt, w0)
        W[i] = np.fft.ifft(xhat * np.conj(daughter))
    return W


def morlet_cwt(
    B: np.ndarray,
    dt: float,
    *,
    w0: float = 6.0,
    dj: float = 1.0 / 12.0,
    s0_factor: float = 2.0,
) -> CWT:
    """Morlet CWT of a (N,3) magnetic field, with COI and scales."""
    B = _fill_gaps(np.asarray(B, dtype=float))
    n = B.shape[0]
    scales, freqs = make_scales(n, dt, dj, s0_factor, w0)
    W = np.stack([cwt_component(B[:, c], dt, scales, w0) for c in range(3)], axis=0)

    # cone of influence (Torrence & Compo 1998, sec. 3g): period = (FF/sqrt2)*d
    t = np.arange(n) * dt
    d = np.minimum(t - t[0], t[-1] - t)
    coi = (fourier_factor(w0) / np.sqrt(2.0)) * d
    coi[coi == 0] = 1e-9
    return CWT(W=W, scales=scales, freqs=freqs, coi=coi, t=t, dt=dt, w0=w0, dj=dj)


# --------------------------------------------------------------------------- #
# Frequency-resolved polarization from the spectral matrix (bonus, slide 9)
# --------------------------------------------------------------------------- #
@dataclass
class Polarization:
    freqs: np.ndarray
    dop: np.ndarray          # degree of polarization (Ns, N)
    ellipticity: np.ndarray  # (Ns, N)


def polarization(cwt: CWT, n_cycles: float = 5.0, wmin: int = 7) -> Polarization:
    r"""Degree of polarization and ellipticity from S_ij = W_i W_j^* (slide 9).

    The spectral matrix is averaged in time with a **scale-dependent** boxcar of
    ``~n_cycles`` wave periods at each scale.  Without this ensemble average the
    instantaneous matrix is rank-1 and the degree of polarization is identically
    1; the running average is what makes DOP < 1 meaningful.

    * **Degree of polarization** ``DOP = sqrt((3 Tr S^2 - (Tr S)^2)/(2 (Tr S)^2))``.
    * **Ellipticity** ``2 |Im S|_off / Tr S`` (0 = linear, 1 = circular); both are
      rotational invariants, so they are frame-independent.
    """
    W = cwt.W                               # (3, Ns, N)
    ns, n = W.shape[1], W.shape[2]
    S = np.einsum("iSN,jSN->ijSN", W, np.conj(W))     # (3,3,Ns,N) complex
    periods = 1.0 / cwt.freqs
    Ssm = np.empty_like(S)
    for s in range(ns):
        w = int(np.clip(round(n_cycles * periods[s] / cwt.dt), wmin, n))
        k = np.ones(w) / w
        for i in range(3):
            for j in range(3):
                Ssm[i, j, s] = (np.convolve(S[i, j, s].real, k, "same")
                                + 1j * np.convolve(S[i, j, s].imag, k, "same"))
    trS = np.einsum("iiSN->SN", Ssm).real
    trS2 = np.einsum("ijSN,jiSN->SN", Ssm, Ssm).real
    with np.errstate(divide="ignore", invalid="ignore"):
        dop = np.sqrt(np.clip((3.0 * trS2 - trS ** 2) / (2.0 * trS ** 2), 0, 1))
        im = (Ssm[0, 1].imag ** 2 + Ssm[1, 2].imag ** 2 + Ssm[2, 0].imag ** 2)
        ellip = 2.0 * np.sqrt(np.clip(im, 0, None)) / trS
    return Polarization(freqs=cwt.freqs, dop=dop, ellipticity=np.clip(ellip, 0, 1))


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _fill_gaps_1d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    bad = ~np.isfinite(x)
    if bad.any():
        idx = np.arange(x.size)
        x = x.copy()
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x


def _fill_gaps(B: np.ndarray) -> np.ndarray:
    return np.column_stack([_fill_gaps_1d(B[:, c]) for c in range(B.shape[1])])


def variance_check(B: np.ndarray, dt: float, fft: FFTSpectrum, cwt: CWT) -> dict:
    """Parseval / energy-conservation tests for both estimators.

    * **FFT (periodogram)** — the full-record one-sided periodogram must
      integrate to the field variance:  ``∫ S(f) df = var(B)``.
    * **CWT (Torrence & Compo eq. 14)** — the exact wavelet reconstruction of
      the variance,  ``σ² = (δj δt / C_δ) Σ_j ⟨|W(s_j)|²⟩ / s_j``.

    (Welch's integral is intentionally *not* compared to the full variance: its
    finite segment length removes power below ``1/seg_s``.)
    """
    B = _fill_gaps(np.asarray(B, dtype=float))
    var_total = float(np.sum(np.var(B, axis=0)))

    # FFT periodogram (full record, mean-removed) -> clean Parseval test
    fr, pr = signal.periodogram(B, fs=1.0 / dt, window="boxcar",
                                detrend="constant", scaling="density", axis=0)
    fft_int = float(np.trapezoid(pr.sum(axis=1), fr))

    # CWT exact reconstruction (Torrence & Compo 1998, eq. 14)
    cdelta = C_DELTA.get(cwt.w0, 0.776)
    mean_pow = cwt.trace_power.mean(axis=1)        # <|W(s_j)|^2> over time, summed comps
    cwt_recon = (cwt.dj * cwt.dt / cdelta) * np.sum(mean_pow / cwt.scales)

    return {
        "var_total": var_total,
        "fft_periodogram_integral": fft_int,
        "cwt_reconstruction": float(cwt_recon),
        "fft_ratio": fft_int / var_total,
        "cwt_ratio": float(cwt_recon) / var_total,
    }
