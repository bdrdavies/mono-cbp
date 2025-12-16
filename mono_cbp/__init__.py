"""
mono-cbp: A pipeline for identifying transits of circumbinary planets in TESS light curves.

This package provides tools for:
- Eclipse masking in eclipsing binary light curves
- Transit detection and characterization
- Model comparison for event classification
- Injection-retrieval analysis for completeness studies
"""

from .eclipse_masking import EclipseMasker
from .transit_finding import TransitFinder
from .model_comparison import ModelComparator
from .injection_retrieval import TransitInjector
from .pipeline import MonoCBPPipeline

__version__ = "0.1.2"
__author__ = "Benjamin Davies"
__all__ = [
    "EclipseMasker",
    "TransitFinder",
    "ModelComparator",
    "TransitInjector",
    "MonoCBPPipeline",
]
