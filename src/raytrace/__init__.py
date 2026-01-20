"""
Real-ray aberration calculation module using Optiland.

This module provides post-ranking real-ray aberration calculation
for cemented and air-spaced doublets using the Optiland ray tracing library.
"""

from .optiland_interface import (
    compute_realray_aberrations,
    RealRayResults,
    OPTILAND_AVAILABLE,
)

__all__ = [
    "compute_realray_aberrations",
    "RealRayResults",
    "OPTILAND_AVAILABLE",
]
