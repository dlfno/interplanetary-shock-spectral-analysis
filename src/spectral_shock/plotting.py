r"""
Figures
=======

Publication-style figures for the README:

* :func:`plot_overview`     — magnetic-field time series + shock + intervals.
* :func:`plot_psd`          — FFT vs CWT power spectral density (up/down).
* :func:`plot_scalograms`   — 2-D CWT scalograms (up/down) + shock crossing.
* :func:`plot_mva`          — MVA hodograms (slide-11 style).
* :func:`plot_polarization` — frequency-resolved DOP/ellipticity (bonus).
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from . import mva as mva_mod
from . import spectral as sp
from .config import EVENT, PARAMS, ShockEvent
from .data import MagData

mpl.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.axisbelow": True,
    "legend.framealpha": 0.9,
    "legend.fontsize": 9,
    "mathtext.default": "regular",
})

C_UP = "#1f6feb"      # upstream colour (blue)
C_DN = "#d62728"      # downstream colour (red)
C_SHOCK = "#111111"
SCAL_CMAP = "turbo"


# --------------------------------------------------------------------------- #
def plot_overview(md: MagData, event: ShockEvent = EVENT, path=None):
    """B_GSE components and |B| vs time, with shock line and analysis windows."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.4]})
    m = md.minutes
    labels = ("$B_x$", "$B_y$", "$B_z$")
    colors = ("#d62728", "#2ca02c", "#1f6feb")
    for c in range(3):
        axes[0].plot(m, md.B[:, c], lw=0.6, color=colors[c], label=labels[c])
    axes[0].set_ylabel("B (GSE)  [nT]")
    axes[0].legend(ncol=3, loc="upper left")

    axes[1].plot(m, md.Bmag, lw=0.7, color="#333333")
    axes[1].set_ylabel("|B|  [nT]")
    axes[1].set_xlabel("Minutes from shock  (%s UT)" % event.t0.strftime("%Y-%m-%d %H:%M:%S"))

    for ax in axes:
        ax.axvline(0, color=C_SHOCK, lw=1.4, ls="-")
        ax.axvspan(*event.upstream_min, color=C_UP, alpha=0.12)
        ax.axvspan(*event.downstream_min, color=C_DN, alpha=0.12)
    # labels
    ymid = axes[0].get_ylim()
    axes[0].text(np.mean(event.upstream_min), ymid[1], " upstream",
                 color=C_UP, va="top", ha="center", fontsize=10, fontweight="bold")
    axes[0].text(np.mean(event.downstream_min), ymid[1], "downstream ",
                 color=C_DN, va="top", ha="center", fontsize=10, fontweight="bold")
    axes[0].text(0, ymid[1], " shock", color=C_SHOCK, va="top", ha="left", fontsize=9)

    src = "Wind/MFI (WI_H2_MFI)" if md.source == "wind" else "SYNTHETIC fallback"
    fig.suptitle(f"{event.spacecraft} fast-forward IP shock — {src}\n"
                 fr"$\theta_{{Bn}}\approx{event.theta_Bn_deg:.0f}°$, "
                 fr"$V_{{sh}}\approx{event.v_shock_kms:.0f}$ km/s, "
                 fr"$M_{{fast}}\approx{event.mach_fast:.2f}$, "
                 fr"$\beta_{{up}}\approx{event.beta_up:.2f}$",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    if path:
        fig.savefig(path)
        plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def _ref_slope(ax, f, anchor_psd, slope, label, x0, x1, color="0.4"):
    """Draw a power-law reference line f^slope."""
    xx = np.array([x0, x1])
    yy = anchor_psd * (xx / x0) ** slope
    ax.plot(xx, yy, ls="--", lw=1.1, color=color, label=label)


def plot_psd(intervals: dict[str, MagData], path=None):
    """FFT (periodogram + Welch) vs CWT global PSD, upstream & downstream."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0), sharey=True)
    order = ["upstream", "downstream"]
    cols = {"upstream": C_UP, "downstream": C_DN}
    for ax, name in zip(axes, order):
        seg = intervals[name]
        B, dt = seg.B, seg.dt
        fft = sp.fft_psd(B, dt, seg_s=PARAMS.welch_seg_s, overlap=PARAMS.welch_overlap,
                         window=PARAMS.window, detrend=PARAMS.detrend)
        cwt = sp.morlet_cwt(B, dt, w0=PARAMS.w0, dj=PARAMS.dj, s0_factor=PARAMS.s0_factor)
        f_cwt, psd_cwt = cwt.global_psd(mask_coi=True)
        o = np.argsort(f_cwt)

        # raw periodogram (grey), Welch (colour), CWT global (green)
        ax.loglog(fft.f_raw[1:], fft.psd_total_raw[1:], color="0.7", lw=0.5,
                  alpha=0.8, label="FFT periodogram (raw)")
        ax.loglog(fft.f[1:], fft.psd_total[1:], color=cols[name], lw=1.8,
                  label="FFT Welch (avg)")
        ax.loglog(f_cwt[o], psd_cwt[o], color="#2ca02c", lw=1.8, label="CWT global")

        # Kolmogorov reference
        fa = 0.05
        idx = np.argmin(np.abs(fft.f - fa))
        _ref_slope(ax, fft.f, fft.psd_total[idx] * 1.5, -5 / 3,
                   r"$f^{-5/3}$", fa, 2.0)

        ax.axvline(seg.fs / 2, color="0.5", ls=":", lw=1)
        ax.text(seg.fs / 2, ax.get_ylim()[1], " Nyquist", rotation=90,
                va="top", ha="right", color="0.4", fontsize=8)
        ax.set_title(f"{name.capitalize()}  ({_window_label(name)})")
        ax.set_xlabel("Frequency  [Hz]")
        ax.legend(loc="lower left")
        ax.grid(which="both", alpha=0.2)
    axes[0].set_ylabel(r"PSD of $\mathbf{B}$  [nT$^2$/Hz]")
    fig.suptitle("Power spectral density — FFT vs CWT", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if path:
        fig.savefig(path)
        plt.close(fig)
    return fig


def _window_label(name: str) -> str:
    rng = EVENT.upstream_min if name == "upstream" else EVENT.downstream_min
    return f"[{rng[0]:+.0f}, {rng[1]:+.0f}] min"


# --------------------------------------------------------------------------- #
def _draw_scalogram(ax, cwt: sp.CWT, t_min, title, cmap=SCAL_CMAP,
                    vmin=None, vmax=None):
    S = cwt.scalogram                                   # (Ns, N) nT^2/Hz
    logS = np.log10(np.clip(S, 1e-12, None))
    pcm = ax.pcolormesh(t_min, cwt.freqs, logS, shading="auto",
                        cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yscale("log")
    # cone of influence
    coi_f = 1.0 / cwt.coi
    ax.plot(t_min, coi_f, color="w", lw=1.2)
    ax.fill_between(t_min, coi_f, cwt.freqs.min(), color="w", alpha=0.35,
                    hatch="xx", edgecolor="0.6", lw=0.0)
    ax.set_ylim(cwt.freqs.min(), cwt.freqs.max())
    ax.set_title(title)
    return pcm


def plot_scalograms(md: MagData, intervals: dict[str, MagData],
                    event: ShockEvent = EVENT, path=None):
    """2-D CWT scalograms: full shock crossing (top) + up/down (bottom)."""
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.32, wspace=0.12)

    # --- full crossing (use full downloaded record clipped to +-11 min) ----
    full = md.slice_minutes(event.upstream_min[0], event.downstream_min[1])
    cwt_full = sp.morlet_cwt(full.B, full.dt, w0=PARAMS.w0, dj=PARAMS.dj,
                             s0_factor=PARAMS.s0_factor)
    axf = fig.add_subplot(gs[0, :])
    # shared colour scale from the full record (robust percentiles)
    Sfull = np.log10(np.clip(cwt_full.scalogram, 1e-12, None))
    vmin, vmax = np.nanpercentile(Sfull, [5, 99.5])
    pcm = _draw_scalogram(axf, cwt_full, full.minutes,
                          "CWT scalogram across the shock crossing",
                          vmin=vmin, vmax=vmax)
    axf.axvline(0, color="k", lw=1.3)
    axf.axvspan(*event.upstream_min, color=C_UP, alpha=0.0)
    axf.set_ylabel("Frequency  [Hz]")
    axf.set_xlabel("Minutes from shock")
    cb = fig.colorbar(pcm, ax=axf, pad=0.01)
    cb.set_label(r"$\log_{10}$ PSD  [nT$^2$/Hz]")

    # --- per-interval -------------------------------------------------------
    for ax_i, name in zip([fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])],
                          ["upstream", "downstream"]):
        seg = intervals[name]
        cwt = sp.morlet_cwt(seg.B, seg.dt, w0=PARAMS.w0, dj=PARAMS.dj,
                            s0_factor=PARAMS.s0_factor)
        pcm = _draw_scalogram(ax_i, cwt, seg.minutes,
                              f"{name.capitalize()}  {_window_label(name)}",
                              vmin=vmin, vmax=vmax)
        ax_i.set_xlabel("Minutes from shock")
        if name == "upstream":
            ax_i.set_ylabel("Frequency  [Hz]")
    fig.colorbar(pcm, ax=fig.axes[1:], pad=0.01).set_label(r"$\log_{10}$ PSD  [nT$^2$/Hz]")
    fig.suptitle("Time–frequency wave activity (Morlet CWT)", fontsize=13)
    if path:
        fig.savefig(path)
        plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def _hodogram(ax, x, y, t_s, xlabel, ylabel, cmap="turbo"):
    pts = np.array([x, y]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=cmap, array=t_s[:-1], lw=0.7, alpha=0.85)
    ax.add_collection(lc)
    lim = 1.08 * max(np.abs(x).max(), np.abs(y).max())
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axhline(0, color="0.7", lw=0.6); ax.axvline(0, color="0.7", lw=0.6)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    return lc


def plot_mva(intervals: dict[str, MagData], event: ShockEvent = EVENT, path=None):
    """MVA hodograms in the variance frame for upstream & downstream."""
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.6))
    rows = {"downstream": 0, "upstream": 1}     # slide-11 order: downstream top
    rcolor = {"upstream": C_UP, "downstream": C_DN}
    last_lc = None
    for name, r in rows.items():
        seg = intervals[name]
        res = mva_mod.mva(seg.B, fs=seg.fs, band_hz=PARAMS.mva_band_hz)
        t = seg.seconds - seg.seconds[0]
        b1, b2, b3 = res.B_mva[:, 0], res.B_mva[:, 1], res.B_mva[:, 2]
        # left: max-min plane (dB1 vs dB3); right: max-int plane (dB1 vs dB2)
        _hodogram(axes[r, 0], b3, b2, t, r"$\Delta B_3$ (min) [nT]", r"$\Delta B_2$ (int) [nT]")
        last_lc = _hodogram(axes[r, 1], b1, b2, t, r"$\Delta B_1$ (max) [nT]", r"$\Delta B_2$ (int) [nT]")
        axes[r, 0].text(0.04, 0.96, fr"$\lambda_2/\lambda_3={res.r23:.2f}$",
                        transform=axes[r, 0].transAxes, va="top", fontsize=11,
                        fontweight="bold", color=rcolor[name])
        axes[r, 1].text(0.04, 0.96,
                        fr"$\lambda_1/\lambda_2={res.r12:.2f}$" + "\n"
                        fr"$\epsilon={res.ellipticity:.2f}$ ({res.polarization_kind})",
                        transform=axes[r, 1].transAxes, va="top", fontsize=11,
                        fontweight="bold", color=rcolor[name])
        axes[r, 0].set_title(f"{name.capitalize()}  {_window_label(name)}",
                             color=rcolor[name], fontweight="bold", loc="left")
        axes[r, 1].set_title(fr"$\theta_{{kB}}={res.theta_kB():.0f}°$, {res.sense()}",
                             loc="right", fontsize=10)
    cb = fig.colorbar(last_lc, ax=axes.ravel().tolist(), pad=0.02, shrink=0.85)
    cb.set_label("seconds into interval")
    fig.suptitle(f"Minimum Variance Analysis — band-pass {PARAMS.mva_band_hz[0]}–"
                 f"{PARAMS.mva_band_hz[1]} Hz", fontsize=13)
    if path:
        fig.savefig(path)
        plt.close(fig)
    return fig


# --------------------------------------------------------------------------- #
def plot_polarization(intervals: dict[str, MagData], path=None):
    """Frequency-resolved degree of polarization & ellipticity (bonus, slide 9)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex="col")
    for col, name in enumerate(["upstream", "downstream"]):
        seg = intervals[name]
        cwt = sp.morlet_cwt(seg.B, seg.dt, w0=PARAMS.w0, dj=PARAMS.dj,
                            s0_factor=PARAMS.s0_factor)
        pol = sp.polarization(cwt, n_cycles=5.0)
        mask = cwt.coi_mask()
        dop = np.where(mask, np.nan, pol.dop)
        ell = np.where(mask, np.nan, pol.ellipticity)
        for row, (data, label, cmap, vlim) in enumerate([
            (dop, "Degree of polarization", "viridis", (0, 1)),
            (ell, "Ellipticity", "plasma", (0, 1)),
        ]):
            ax = axes[row, col]
            pcm = ax.pcolormesh(seg.minutes, cwt.freqs, data, shading="auto",
                                cmap=cmap, vmin=vlim[0], vmax=vlim[1])
            ax.plot(seg.minutes, 1 / cwt.coi, color="w", lw=1)
            ax.set_yscale("log")
            ax.set_ylim(cwt.freqs.min(), cwt.freqs.max())
            if col == 0:
                ax.set_ylabel("Frequency [Hz]")
            if row == 0:
                ax.set_title(name.capitalize())
            if row == 1:
                ax.set_xlabel("Minutes from shock")
            fig.colorbar(pcm, ax=ax, pad=0.01).set_label(label)
    fig.suptitle("Frequency-resolved polarization (CWT spectral matrix)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if path:
        fig.savefig(path)
        plt.close(fig)
    return fig
