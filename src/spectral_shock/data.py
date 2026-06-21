"""
Data access
===========

Download, cache and load the high-resolution *Wind*/MFI magnetic field for the
configured shock event.

Primary source
--------------
NASA/GSFC **CDAWeb CDAS REST API**.  We request a *subset* of the dataset
``WI_H2_MFI`` (variable ``BGSE``) for the analysis window only, so the download
is a few hundred kB instead of a full-day file::

    https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/
        WI_H2_MFI/data/<start>,<stop>/BGSE?format=cdf

The returned CDF is read with :mod:`cdflib`.  Results are cached as a small
``.npz`` so the analysis reproduces **offline** once the data have been fetched.

Fallback
--------
If the network / CDAWeb is unavailable, :func:`load_event` can synthesise a
*physically flavoured* shock (clearly labelled ``source='synthetic'``) so the
pipeline always runs.  Real data are strongly preferred and used by default.
"""

from __future__ import annotations

import json
import ssl
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

from .config import DATA_DIR, EVENT, ShockEvent, fill_to_nan

CDAS_BASE = "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets"


# --------------------------------------------------------------------------- #
# Container
# --------------------------------------------------------------------------- #
@dataclass
class MagData:
    """A magnetic-field time series in a Cartesian frame (default GSE)."""

    time: np.ndarray          # np.datetime64[ns], shape (N,)
    B: np.ndarray             # shape (N, 3) [nT]
    t0: datetime              # reference epoch (shock time)
    source: str = "wind"      # 'wind' or 'synthetic'
    components: tuple[str, str, str] = ("Bx", "By", "Bz")
    frame: str = "GSE"

    # --- derived ----------------------------------------------------------
    @property
    def seconds(self) -> np.ndarray:
        """Seconds relative to ``t0`` (float)."""
        t0 = np.datetime64(self.t0.replace(tzinfo=None), "ns")
        return (self.time - t0) / np.timedelta64(1, "s")

    @property
    def minutes(self) -> np.ndarray:
        return self.seconds / 60.0

    @property
    def Bmag(self) -> np.ndarray:
        return np.linalg.norm(self.B, axis=1)

    @property
    def dt(self) -> float:
        """Median sampling interval [s]."""
        d = np.diff(self.time) / np.timedelta64(1, "s")
        return float(np.median(d))

    @property
    def fs(self) -> float:
        return 1.0 / self.dt

    def __len__(self) -> int:
        return self.B.shape[0]

    def slice_minutes(self, m0: float, m1: float) -> "MagData":
        mins = self.minutes
        sel = (mins >= m0) & (mins < m1)
        return MagData(
            time=self.time[sel],
            B=self.B[sel],
            t0=self.t0,
            source=self.source,
            components=self.components,
            frame=self.frame,
        )


# --------------------------------------------------------------------------- #
# CDAWeb fetch
# --------------------------------------------------------------------------- #
def _iso_compact(dt: datetime) -> str:
    return dt.strftime("%Y%m%dT%H%M%SZ")


def _cdas_url(event: ShockEvent) -> str:
    start, stop = event.download_window()
    return (
        f"{CDAS_BASE}/{event.dataset}/data/"
        f"{_iso_compact(start)},{_iso_compact(stop)}/{event.variable}?format=cdf"
    )


def _download_cdf(event: ShockEvent, dest: Path, timeout: int = 120) -> Path:
    """Resolve the CDAS REST request to a CDF file and download it."""
    ctx = ssl.create_default_context()
    url = _cdas_url(event)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        meta = json.load(r)
    files = meta.get("FileDescription", [])
    if not files:
        raise RuntimeError(
            f"CDAS returned no files. Status={meta.get('Status')} "
            f"Error={meta.get('Error')} Warning={meta.get('Warning')}"
        )
    cdf_url = files[0]["Name"]
    with urllib.request.urlopen(cdf_url, timeout=timeout, context=ctx) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def _read_cdf(path: Path, event: ShockEvent) -> MagData:
    import cdflib

    cdf = cdflib.CDF(str(path))
    epoch = cdf.varget(event.epoch_var)
    bgse = cdf.varget(event.variable)
    time = cdflib.cdfepoch.to_datetime(epoch).astype("datetime64[ns]")
    B = fill_to_nan(np.asarray(bgse, dtype=float))
    return MagData(time=time, B=B, t0=event.t0, source="wind")


# --------------------------------------------------------------------------- #
# Caching
# --------------------------------------------------------------------------- #
def _npz_path(event: ShockEvent) -> Path:
    return DATA_DIR / f"{event.spacecraft}_{event.dataset}_{event.slug}.npz"


def _cdf_path(event: ShockEvent) -> Path:
    return DATA_DIR / f"{event.spacecraft}_{event.dataset}_{event.slug}.cdf"


def _save_npz(path: Path, md: MagData) -> None:
    np.savez_compressed(
        path,
        time=md.time.astype("datetime64[ns]").astype("int64"),
        B=md.B,
        t0=np.datetime64(md.t0.replace(tzinfo=None), "ns").astype("int64"),
        source=md.source,
        frame=md.frame,
        components=np.array(md.components),
    )


def _load_npz(path: Path, event: ShockEvent) -> MagData:
    z = np.load(path, allow_pickle=False)
    return MagData(
        time=z["time"].astype("datetime64[ns]"),
        B=z["B"],
        t0=event.t0,
        source=str(z["source"]),
        components=tuple(z["components"]),
        frame=str(z["frame"]),
    )


# --------------------------------------------------------------------------- #
# Synthetic fallback
# --------------------------------------------------------------------------- #
def synthesize(event: ShockEvent, fs: float = 10.8696, seed: int = 20220819) -> MagData:
    """A physically flavoured synthetic shock (clearly labelled).

    Upstream: quasi-parallel foreshock — quiet mean field with a transverse,
    nearly *circularly polarised* ULF wave (~0.05 Hz) plus turbulence.
    Downstream: compressed (~2x) mean field with larger-amplitude, more
    *elliptical/compressive* turbulence.  Used only when CDAWeb is unreachable.
    """
    rng = np.random.default_rng(seed)
    start, stop = event.download_window()
    n = int((stop - start).total_seconds() * fs)
    t_s = np.arange(n) / fs
    time = (np.datetime64(start.replace(tzinfo=None), "ns")
            + (t_s * 1e9).astype("timedelta64[ns]"))
    t_rel = t_s - (event.t0 - start).total_seconds()      # seconds from shock
    ramp = 0.5 * (1.0 + np.tanh(t_rel / 3.0))             # smooth shock ramp

    # mean field rotates and compresses across the shock
    b_up = np.array([3.0, -3.2, 2.0])
    b_dn = np.array([5.0, -8.0, 6.0])
    mean = b_up[None, :] * (1 - ramp)[:, None] + b_dn[None, :] * ramp[:, None]

    # circularly polarised transverse wave upstream (in a plane ~perp to B)
    f_ulf = 0.05
    e1 = np.array([0.0, 0.6, 0.8]); e1 /= np.linalg.norm(e1)
    e2 = np.array([1.0, 0.0, 0.0])
    a_up = 1.1 * (1 - ramp)
    wave_up = (a_up[:, None] * (np.cos(2 * np.pi * f_ulf * t_s)[:, None] * e1
                                + np.sin(2 * np.pi * f_ulf * t_s)[:, None] * e2))
    # elliptical, compressive wave downstream (mostly along one axis)
    f_dn = 0.12
    a_dn = 2.4 * ramp
    wave_dn = a_dn[:, None] * np.column_stack([
        1.0 * np.cos(2 * np.pi * f_dn * t_s),
        0.35 * np.sin(2 * np.pi * f_dn * t_s),
        0.2 * np.cos(2 * np.pi * f_dn * t_s + 0.7),
    ])

    # turbulence: stronger downstream, with a ~ f^-1.6 colouring
    noise = rng.standard_normal((n, 3))
    color = np.fft.rfft(noise, axis=0)
    freqs = np.fft.rfftfreq(n, d=1 / fs)
    freqs[0] = freqs[1]
    color *= (freqs[:, None] ** (-0.8))
    turb = np.fft.irfft(color, n=n, axis=0)
    turb /= turb.std(axis=0, keepdims=True)
    turb *= (0.35 + 1.6 * ramp)[:, None]

    B = mean + wave_up + wave_dn + turb
    return MagData(time=time, B=B, t0=event.t0, source="synthetic")


# --------------------------------------------------------------------------- #
# Public loader
# --------------------------------------------------------------------------- #
def load_event(
    event: ShockEvent = EVENT,
    *,
    use_cache: bool = True,
    allow_download: bool = True,
    allow_synthetic: bool = True,
    force_synthetic: bool = False,
    verbose: bool = True,
) -> MagData:
    """Load the magnetic field for ``event`` (cache -> download -> synthetic).

    Synthetic data are *never* written to the real-data cache file, so forcing
    or falling back to the synthetic model cannot poison the committed cache.
    """
    if force_synthetic:
        if verbose:
            print("[data] >>> forcing SYNTHETIC data (no download) <<<")
        return synthesize(event)

    npz = _npz_path(event)
    if use_cache and npz.exists():
        if verbose:
            print(f"[data] loading cached {npz.name}")
        return _load_npz(npz, event)

    if allow_download:
        try:
            cdf = _cdf_path(event)
            if verbose:
                print(f"[data] downloading {event.dataset} {event.variable} from CDAWeb ...")
            _download_cdf(event, cdf)
            md = _read_cdf(cdf, event)
            _save_npz(npz, md)            # only *real* data are cached
            if verbose:
                print(f"[data] fetched {len(md)} samples @ {md.fs:.3f} Hz "
                      f"({md.time[0]} .. {md.time[-1]})")
            return md
        except Exception as exc:  # noqa: BLE001
            if verbose:
                print(f"[data] download failed ({exc!r})")
            if not allow_synthetic:
                raise

    if allow_synthetic:
        if verbose:
            print("[data] >>> using SYNTHETIC fallback data (not cached) <<<")
        return synthesize(event)

    raise RuntimeError("No data available (cache miss, download disabled).")


def get_intervals(md: MagData, event: ShockEvent = EVENT) -> dict[str, MagData]:
    """Return {'upstream':..., 'downstream':...} sliced views."""
    up = md.slice_minutes(*event.upstream_min)
    dn = md.slice_minutes(*event.downstream_min)
    return {"upstream": up, "downstream": dn}


if __name__ == "__main__":
    md = load_event()
    print(f"source={md.source} N={len(md)} fs={md.fs:.4f} Hz dt={md.dt:.4f} s")
    print(f"|B| upstream window vs downstream window:")
    iv = get_intervals(md)
    for k, v in iv.items():
        print(f"  {k:<11s} N={len(v):5d}  <|B|>={np.nanmean(v.Bmag):6.2f} nT  "
              f"std|B|={np.nanstd(v.Bmag):5.2f}")
