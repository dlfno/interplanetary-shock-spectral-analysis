"""
Offline unit tests for the spectral / MVA core (no network required).

Run with:  python -m pytest -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from spectral_shock import mva, spectral  # noqa: E402


def _synth(n=8192, dt=0.1, freqs=(0.05, 0.3), amps=(2.0, 1.0), seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * dt
    B = rng.standard_normal((n, 3)) * 0.3
    for f, a in zip(freqs, amps):
        B[:, 0] += a * np.cos(2 * np.pi * f * t)
        B[:, 1] += a * np.sin(2 * np.pi * f * t)   # circular in x-y
    return B, dt


def test_fft_parseval():
    """Full-record periodogram integrates to the field variance (Parseval)."""
    B, dt = _synth()
    fft = spectral.fft_psd(B, dt)
    cwt = spectral.morlet_cwt(B, dt)
    chk = spectral.variance_check(B, dt, fft, cwt)
    assert abs(chk["fft_ratio"] - 1.0) < 1e-2


def test_cwt_reconstruction():
    """Torrence & Compo wavelet reconstruction recovers the variance to ~10%."""
    B, dt = _synth()
    fft = spectral.fft_psd(B, dt)
    cwt = spectral.morlet_cwt(B, dt)
    chk = spectral.variance_check(B, dt, fft, cwt)
    assert 0.85 < chk["cwt_ratio"] < 1.10


def test_fourier_factor():
    """Morlet Fourier factor for omega0 = 6 is ~1.033 (scale ~ period)."""
    assert abs(spectral.fourier_factor(6.0) - 1.0330) < 1e-3


def test_cwt_psd_matches_fft():
    """The CWT global spectrum and FFT Welch agree at an injected spectral peak."""
    B, dt = _synth(freqs=(0.1,), amps=(3.0,))
    fft = spectral.fft_psd(B, dt, seg_s=200.0)
    cwt = spectral.morlet_cwt(B, dt)
    f_c, psd_c = cwt.global_psd(mask_coi=True)
    fpk_fft = fft.f[1:][np.argmax(fft.psd_total[1:])]
    fpk_cwt = f_c[np.nanargmax(psd_c)]
    assert abs(fpk_fft - 0.1) < 0.02
    assert abs(fpk_cwt - 0.1) < 0.03


def test_mva_circular():
    """A circular x-y wave -> lambda1/lambda2 ~ 1 and a small lambda3."""
    B, dt = _synth(freqs=(0.05,), amps=(3.0,), seed=1)
    res = mva.mva(B, fs=1 / dt, band_hz=(0.02, 0.2))
    assert res.r12 < 1.5                  # near-circular
    assert res.polarization_kind == "circular"
    # minimum-variance direction ~ z-axis (plane is x-y)
    assert abs(abs(res.normal[2]) - 1.0) < 0.15


def test_mva_linear():
    """A linear wave (single axis) -> large lambda1/lambda2."""
    n, dt = 8192, 0.1
    t = np.arange(n) * dt
    B = np.zeros((n, 3))
    B[:, 0] = 3.0 * np.cos(2 * np.pi * 0.05 * t)
    B += np.random.default_rng(2).standard_normal((n, 3)) * 0.05
    res = mva.mva(B, fs=1 / dt, band_hz=(0.02, 0.2))
    assert res.r12 > 4.0
    assert res.polarization_kind == "linear"
