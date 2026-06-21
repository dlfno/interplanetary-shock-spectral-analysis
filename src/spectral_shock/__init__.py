"""
spectral_shock
==============

Spectral analysis of an interplanetary (collisionless) shock observed by the
NASA *Wind* spacecraft, for the UNAM course *Space Science and Data Analysis
(3. Spectral Analysis)*.

Modules
-------
config    : event definition, intervals and physical constants.
data      : download / cache / load Wind MFI high-resolution magnetic field.
spectral  : FFT power spectral density and a Torrence & Compo (1998) Morlet CWT.
mva       : minimum variance analysis and wave-polarization diagnostics.
plotting  : publication-style figures used in the README.
"""

from . import config, data, spectral, mva  # noqa: F401

__all__ = ["config", "data", "spectral", "mva"]
__version__ = "1.0.0"
