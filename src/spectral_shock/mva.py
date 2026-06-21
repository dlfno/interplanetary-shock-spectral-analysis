r"""
Minimum Variance Analysis (MVA)
===============================

Wave-property analysis following Sonnerup & Scheible (in Paschmann & Daly,
*Analysis Methods for Multi-Spacecraft Data*, ISSI 1998).

Given a magnetic-field interval :math:`\mathbf{B}^{(m)}`, build the variance
matrix

.. math::  M_{ij} = \langle B_i B_j\rangle - \langle B_i\rangle\langle B_j\rangle ,

and diagonalise it.  The eigenvalues :math:`\lambda_1\ge\lambda_2\ge\lambda_3`
(maximum, intermediate, minimum) with eigenvectors
:math:`\hat{\mathbf{x}}_1,\hat{\mathbf{x}}_2,\hat{\mathbf{x}}_3` define the
variance ellipsoid.  :math:`\hat{\mathbf{x}}_3` is the minimum-variance
direction — the estimated wave-normal / propagation direction.

Interpretation of the eigenvalue ratios
---------------------------------------
* :math:`\lambda_2/\lambda_3 \gg 1`  → the minimum-variance direction
  (polarization plane / wave normal) is **well defined**.
* :math:`\lambda_1/\lambda_2 \approx 1` → **circular** polarization;
  :math:`\lambda_1/\lambda_2 \gg 1` → **linear**; in between → **elliptical**.
* A common reliability threshold is :math:`\lambda_2/\lambda_3 \gtrsim 5\text{–}10`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal


@dataclass
class MVAResult:
    eigvals: np.ndarray        # (3,) lambda1>=lambda2>=lambda3
    eigvecs: np.ndarray        # (3,3) columns = x1 (max), x2 (int), x3 (min)
    B_mva: np.ndarray          # (N,3) fluctuation field in MVA frame (dB1,dB2,dB3)
    B_mean: np.ndarray         # (3,) mean field in input frame
    band_hz: tuple[float, float] | None
    n_used: int

    # --- ratios -----------------------------------------------------------
    @property
    def lam1(self) -> float: return float(self.eigvals[0])
    @property
    def lam2(self) -> float: return float(self.eigvals[1])
    @property
    def lam3(self) -> float: return float(self.eigvals[2])
    @property
    def r12(self) -> float: return self.lam1 / self.lam2     # linear vs circular
    @property
    def r23(self) -> float: return self.lam2 / self.lam3     # quality of normal

    @property
    def normal(self) -> np.ndarray:
        """Minimum-variance direction (wave normal), in the input frame."""
        return self.eigvecs[:, 2]

    def theta_kB(self) -> float:
        """Angle [deg] between the wave normal and the mean field."""
        n = self.normal
        b = self.B_mean / np.linalg.norm(self.B_mean)
        c = abs(float(np.dot(n, b)))
        return float(np.degrees(np.arccos(np.clip(c, -1, 1))))

    @property
    def ellipticity(self) -> float:
        r"""Ellipticity of the hodogram in the maximum-variance plane.

        Estimated from the variance ratio, :math:`\sqrt{\lambda_2/\lambda_1}`:
        1 = circular, 0 = linear.
        """
        return float(np.sqrt(self.lam2 / self.lam1))

    def sense(self) -> str:
        """Sense of rotation of dB in the (x1,x2) plane w.r.t. the mean field.

        Right-handed / left-handed about +B is inferred from the sign of the
        time-averaged z-component of dB1 x dB2 projected on the normal.
        """
        b1 = self.B_mva[:, 0]
        b2 = self.B_mva[:, 1]
        # signed area swept (z of cross product of successive vectors)
        cross = b1[:-1] * b2[1:] - b2[:-1] * b1[1:]
        s = np.sum(cross)
        # relate rotation sense to field direction (x3 vs +B)
        hand = np.sign(np.dot(self.normal, self.B_mean))
        val = s * hand
        if abs(s) < 1e-6:
            return "indeterminate"
        return "right-handed" if val < 0 else "left-handed"

    @property
    def polarization_kind(self) -> str:
        if self.r12 < 1.5:
            return "circular"
        if self.r12 > 4.0:
            return "linear"
        return "elliptical"


def bandpass(B: np.ndarray, fs: float, band: tuple[float, float], order: int = 4) -> np.ndarray:
    """Zero-phase Butterworth band-pass of a (N,3) field."""
    lo, hi = band
    nyq = fs / 2.0
    hi = min(hi, 0.99 * nyq)
    sos = signal.butter(order, [lo / nyq, hi / nyq], btype="band", output="sos")
    return np.column_stack([signal.sosfiltfilt(sos, B[:, c]) for c in range(3)])


def mva(
    B: np.ndarray,
    fs: float | None = None,
    band_hz: tuple[float, float] | None = None,
) -> MVAResult:
    """Minimum variance analysis of a (N,3) magnetic field.

    If ``band_hz`` is given the field is band-pass filtered first, so the
    analysis targets the wave activity in that band.
    """
    B = np.asarray(B, dtype=float)
    B = _fill_gaps(B)
    B_mean = B.mean(axis=0)

    work = B
    if band_hz is not None:
        if fs is None:
            raise ValueError("fs required for band-pass filtering")
        work = bandpass(B, fs, band_hz)

    dB = work - work.mean(axis=0)
    # variance (covariance) matrix M_ij = <dB_i dB_j>
    M = (dB.T @ dB) / dB.shape[0]
    eigvals, eigvecs = np.linalg.eigh(M)          # ascending
    order = np.argsort(eigvals)[::-1]             # -> descending (max, int, min)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # enforce a right-handed eigenvector triad
    if np.dot(np.cross(eigvecs[:, 0], eigvecs[:, 1]), eigvecs[:, 2]) < 0:
        eigvecs[:, 2] *= -1.0
    # point the minimum-variance direction along +B (convention)
    if np.dot(eigvecs[:, 2], B_mean) < 0:
        eigvecs[:, 1] *= -1.0
        eigvecs[:, 2] *= -1.0

    B_mva = dB @ eigvecs                          # (N,3) in (x1,x2,x3)
    return MVAResult(eigvals=eigvals, eigvecs=eigvecs, B_mva=B_mva,
                     B_mean=B_mean, band_hz=band_hz, n_used=B.shape[0])


def _fill_gaps(B: np.ndarray) -> np.ndarray:
    out = np.array(B, dtype=float)
    for c in range(out.shape[1]):
        x = out[:, c]
        bad = ~np.isfinite(x)
        if bad.any():
            idx = np.arange(x.size)
            x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
            out[:, c] = x
    return out


def summary(res: MVAResult) -> str:
    return (
        f"lambda = [{res.lam1:.3g}, {res.lam2:.3g}, {res.lam3:.3g}] nT^2 | "
        f"l1/l2={res.r12:.2f} l2/l3={res.r23:.2f} | "
        f"ellip={res.ellipticity:.2f} ({res.polarization_kind}, {res.sense()}) | "
        f"theta_kB={res.theta_kB():.1f} deg | "
        f"n_min={np.array2string(res.normal, precision=3)}"
    )
