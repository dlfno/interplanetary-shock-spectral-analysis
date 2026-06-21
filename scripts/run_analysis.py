#!/usr/bin/env python3
"""
Run the full spectral-analysis pipeline for the configured shock event.

Steps
-----
1. Load Wind/MFI high-resolution magnetic field (download + cache).
2. Time-series overview figure.
3. FFT vs CWT power spectral density (upstream & downstream).
4. 2-D CWT scalograms (shock crossing + each interval).
5. Minimum variance analysis hodograms.
6. Frequency-resolved polarization (bonus).
7. Write results table to ``figures/results.md`` and ``results.json``.

Usage
-----
    python scripts/run_analysis.py [--synthetic] [--no-cache]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from spectral_shock import data, mva, plotting, spectral          # noqa: E402
from spectral_shock.config import EVENT, FIG_DIR, PARAMS          # noqa: E402


def fmt(x, n=3):
    return f"{x:.{n}g}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--synthetic", action="store_true",
                    help="force the synthetic fallback (no download)")
    ap.add_argument("--no-cache", action="store_true", help="ignore cached data")
    args = ap.parse_args()

    print("=" * 72)
    print(f"  {EVENT.spacecraft} interplanetary shock — {EVENT.t0_iso}")
    print(f"  dataset={EVENT.dataset}  var={EVENT.variable}")
    print("=" * 72)

    md = data.load_event(
        use_cache=not args.no_cache,
        allow_download=not args.synthetic,
        allow_synthetic=True,
        force_synthetic=args.synthetic,
    )
    print(f"[load] source={md.source}  N={len(md)}  fs={md.fs:.4f} Hz  dt={md.dt:.4f} s")

    intervals = data.get_intervals(md)

    # --- figures -----------------------------------------------------------
    print("[fig ] overview ...")
    plotting.plot_overview(md, EVENT, FIG_DIR / "fig01_overview.png")
    print("[fig ] PSD (FFT vs CWT) ...")
    plotting.plot_psd(intervals, FIG_DIR / "fig02_psd_fft_cwt.png")
    print("[fig ] CWT scalograms ...")
    plotting.plot_scalograms(md, intervals, EVENT, FIG_DIR / "fig03_cwt_scalogram.png")
    print("[fig ] MVA hodograms ...")
    plotting.plot_mva(intervals, EVENT, FIG_DIR / "fig04_mva_hodograms.png")
    print("[fig ] polarization ...")
    plotting.plot_polarization(intervals, FIG_DIR / "fig05_polarization.png")

    # --- numbers -----------------------------------------------------------
    results = {
        "event": {
            "spacecraft": EVENT.spacecraft, "dataset": EVENT.dataset,
            "t0": EVENT.t0_iso, "theta_Bn_deg": EVENT.theta_Bn_deg,
            "v_shock_kms": EVENT.v_shock_kms, "mach_fast": EVENT.mach_fast,
            "beta_up": EVENT.beta_up, "catalog_url": EVENT.catalog_url,
            "source": md.source, "fs_hz": md.fs, "n_samples": len(md),
        },
        "intervals": {},
    }
    rows = []
    for name, seg in intervals.items():
        fft = spectral.fft_psd(seg.B, seg.dt, seg_s=PARAMS.welch_seg_s)
        cwt = spectral.morlet_cwt(seg.B, seg.dt, w0=PARAMS.w0, dj=PARAMS.dj,
                                  s0_factor=PARAMS.s0_factor)
        chk = spectral.variance_check(seg.B, seg.dt, fft, cwt)
        res = mva.mva(seg.B, fs=seg.fs, band_hz=PARAMS.mva_band_hz)
        fpk = float(fft.f[1:][np.argmax(fft.psd_total[1:])])
        results["intervals"][name] = {
            "window_min": list(EVENT.upstream_min if name == "upstream"
                               else EVENT.downstream_min),
            "n": len(seg),
            "mean_Bmag_nT": float(np.nanmean(seg.Bmag)),
            "var_total_nT2": chk["var_total"],
            "fft_parseval_ratio": chk["fft_ratio"],
            "cwt_recon_ratio": chk["cwt_ratio"],
            "psd_peak_hz": fpk,
            "mva": {
                "lambda": [res.lam1, res.lam2, res.lam3],
                "l1_l2": res.r12, "l2_l3": res.r23,
                "ellipticity": res.ellipticity,
                "polarization": res.polarization_kind,
                "sense": res.sense(),
                "theta_kB_deg": res.theta_kB(),
                "normal": res.normal.tolist(),
            },
        }
        rows.append((name, seg, chk, res, fpk))

    # console table
    print("\n" + "-" * 72)
    print(f"{'interval':<11s}{'<|B|>':>8s}{'var':>9s}{'l1/l2':>8s}"
          f"{'l2/l3':>8s}{'ellip':>7s}{'pol':>12s}{'th_kB':>7s}")
    print("-" * 72)
    for name, seg, chk, res, fpk in rows:
        print(f"{name:<11s}{np.nanmean(seg.Bmag):>8.2f}{chk['var_total']:>9.2f}"
              f"{res.r12:>8.2f}{res.r23:>8.2f}{res.ellipticity:>7.2f}"
              f"{res.polarization_kind:>12s}{res.theta_kB():>6.0f}°")
    print("-" * 72)

    # --- write results -----------------------------------------------------
    (ROOT / "results.json").write_text(json.dumps(results, indent=2))
    _write_results_md(results, ROOT / "figures" / "results.md")
    print(f"\n[done] figures -> {FIG_DIR}")
    print(f"[done] results -> results.json , figures/results.md")
    return 0


def _write_results_md(results: dict, path: Path):
    e = results["event"]
    lines = [
        "# Results / Resultados",
        "",
        f"**Event:** {e['spacecraft']} interplanetary shock — `{e['t0']}` "
        f"(source: `{e['source']}`, fs = {e['fs_hz']:.3f} Hz, N = {e['n_samples']}).",
        "",
        "| Interval | window | ⟨\\|B\\|⟩ [nT] | var [nT²] | λ₁/λ₂ | λ₂/λ₃ | "
        "ellipticity | polarization | sense | θ_kB | FFT Parseval | CWT recon |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name, iv in results["intervals"].items():
        m = iv["mva"]
        w = iv["window_min"]
        lines.append(
            f"| **{name}** | [{w[0]:+.0f}, {w[1]:+.0f}] min | "
            f"{iv['mean_Bmag_nT']:.2f} | {iv['var_total_nT2']:.2f} | "
            f"{m['l1_l2']:.2f} | {m['l2_l3']:.2f} | {m['ellipticity']:.2f} | "
            f"{m['polarization']} | {m['sense']} | {m['theta_kB_deg']:.0f}° | "
            f"{iv['fft_parseval_ratio']:.3f} | {iv['cwt_recon_ratio']:.3f} |"
        )
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
