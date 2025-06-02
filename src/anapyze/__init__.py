"""
Anapyze: a collection of neuroimaging biomarker tools.

This package exposes:
  - core         (low-level processing utilities)
  - analysis     (statistical routines: two-sample t-tests, correlations)
  - io           (I/O helpers for SPM, CAT12, ADNI, etc.)
"""

from .core import *
from .analysis import *
from .io import *

__all__ = ['core', 'analysis', 'io']