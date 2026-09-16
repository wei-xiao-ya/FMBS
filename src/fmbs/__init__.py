"""FMBS: fused multi-level bioactivity similarity target prediction."""

from .constants import LAYERS
from .training import FMBSParameters, fit_parameters

__all__ = ["FMBSParameters", "LAYERS", "fit_parameters"]
__version__ = "0.1.0"
