"""Utility functions for mono-cbp."""

from .data import (
    load_catalogue,
    get_row,
    bin_to_long_cadence,
    catalogue_to_lc_files,
)
from .eclipses import (
    time_to_phase,
    get_eclipse_mask,
    get_eclipse_indices,
)
from .monofind import (
    monofind,
    get_var_mad,
    split_tol,
)
from .detrending import (
    detrend,
    cosine_detrend,
    slider_detrend,
    run_multi_biweight,
)
from .common import (
    setup_logging,
    get_snr,
)
from .transit_models import (
    create_transit_models,
    save_transit_models,
    load_transit_models,
)

__all__ = [
    "load_catalogue",
    "get_row",
    "bin_to_long_cadence",
    "catalogue_to_lc_files",
    "time_to_phase",
    "get_eclipse_mask",
    "get_eclipse_indices",
    "monofind",
    "get_var_mad",
    "split_tol",
    "detrend",
    "cosine_detrend",
    "slider_detrend",
    "run_multi_biweight",
    "setup_logging",
    "get_snr",
    "create_transit_models",
    "save_transit_models",
    "load_transit_models"
]
