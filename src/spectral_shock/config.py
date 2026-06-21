"""
Event configuration
====================

Single source of truth for the chosen shock event, the analysis intervals and
the data product.  Edit ``EVENT`` to analyse a different shock.

Chosen event
------------
NASA *Wind* spacecraft, fast-forward interplanetary shock on **2022-08-19
16:51:00 UT** (CfA Interplanetary Shock Database #00802).  The shock parameters
quoted below come from the CfA database (Rankine--Hugoniot ``RH08`` solution).

The high-resolution magnetic-field product ``WI_H2_MFI`` (~11 vectors s^-1,
cadence ~0.092 s) is exactly the data set referenced in the lecture notes
("1-hour Wind MFI magnetic field has approximately 39130 data points").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
DATA_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Event definition
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ShockEvent:
    """Definition of a shock event and its up-/downstream analysis windows."""

    name: str
    spacecraft: str
    dataset: str                 # CDAWEB dataset id
    variable: str                # vector magnetic-field variable in the CDF
    epoch_var: str               # time variable in the CDF
    t0: datetime                 # shock arrival time (UT)
    # analysis windows, in minutes relative to t0  (start, stop)
    upstream_min: tuple[float, float]
    downstream_min: tuple[float, float]
    # padding (minutes) downloaded on each side of [upstream_start, downstream_stop]
    pad_min: float = 5.0

    # --- published shock parameters (CfA RH08), for reporting only ----------
    shock_type: str = "FF"       # fast-forward
    theta_Bn_deg: float = 31.3
    theta_Bn_err: float = 11.6
    v_shock_kms: float = 683.4
    mach_fast: float = 1.81
    mach_ms: float = 3.26        # "slow"/magnetosonic Mach from CfA table
    beta_up: float = 0.553
    beta_down: float = 0.796
    n_up_cc: float = 3.1
    v_sw_kms: float = 536.0
    b_up_nT: float = 4.9
    normal_gse: tuple[float, float, float] = (-0.994, 0.114, 0.006)
    catalog_url: str = "https://lweb.cfa.harvard.edu/shocks/wi_data/00802/wi_00802.html"

    # ----- derived helpers --------------------------------------------------
    @property
    def t0_iso(self) -> str:
        return self.t0.strftime("%Y-%m-%dT%H:%M:%SZ")

    def window(self, which: str) -> tuple[datetime, datetime]:
        """Absolute (start, stop) datetimes for 'upstream' or 'downstream'."""
        rng = self.upstream_min if which == "upstream" else self.downstream_min
        return (self._shift(rng[0]), self._shift(rng[1]))

    def download_window(self) -> tuple[datetime, datetime]:
        """Absolute (start, stop) of the full interval to download (with pad)."""
        start = self._shift(self.upstream_min[0] - self.pad_min)
        stop = self._shift(self.downstream_min[1] + self.pad_min)
        return start, stop

    def _shift(self, minutes: float) -> datetime:
        from datetime import timedelta

        return self.t0 + timedelta(minutes=minutes)

    @property
    def slug(self) -> str:
        return self.t0.strftime("%Y%m%d_%H%M%S")


EVENT = ShockEvent(
    name="Wind fast-forward interplanetary shock",
    spacecraft="Wind",
    dataset="WI_H2_MFI",
    variable="BGSE",
    epoch_var="Epoch",
    t0=datetime(2022, 8, 19, 16, 51, 0, tzinfo=timezone.utc),
    upstream_min=(-11.0, -1.0),     # [-11, -1] min  -> 16:40-16:50 UT
    downstream_min=(1.0, 11.0),     # [ +1, +11] min -> 16:52-17:02 UT
    pad_min=6.0,
)


# --------------------------------------------------------------------------- #
# Method parameters
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SpectralParams:
    """Knobs for the spectral analysis."""

    # --- Morlet CWT (Torrence & Compo 1998) --------------------------------
    w0: float = 6.0          # non-dimensional frequency omega_0
    dj: float = 1.0 / 12.0   # scale spacing (sub-octaves); finer than lecture's 1/8
    s0_factor: float = 2.0   # smallest scale s0 = s0_factor * dt
    # --- FFT / Welch -------------------------------------------------------
    welch_seg_s: float = 100.0   # Welch segment length [s]
    welch_overlap: float = 0.5
    window: str = "hann"
    detrend: str = "linear"
    # --- MVA ---------------------------------------------------------------
    # band-pass band (Hz) used to isolate the coherent wave activity before MVA.
    # 0.02-0.25 Hz (periods ~4-50 s) targets the ULF foreshock / whistler band.
    mva_band_hz: tuple[float, float] = (0.02, 0.25)


PARAMS = SpectralParams()


# Physical constants (SI unless noted)
MU0 = 4.0e-7 * np.pi          # vacuum permeability [H/m]
M_PROTON = 1.67262192e-27     # proton mass [kg]
KB = 1.380649e-23             # Boltzmann constant [J/K]


def fill_to_nan(arr: np.ndarray, fill_threshold: float = -1.0e30) -> np.ndarray:
    """Replace CDF fill values with NaN."""
    out = np.array(arr, dtype=float)
    out[out <= fill_threshold] = np.nan
    return out
