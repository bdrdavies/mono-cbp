"""Shared utilities for parsing light curve filenames and loading light curve files.

These functions are the single implementation used by EclipseMasker, TransitFinder
and TransitInjector, which previously each carried their own (subtly different) copies.
"""

import os
import logging
import numpy as np

logger = logging.getLogger('mono_cbp.utils.io')

DEFAULT_NPZ_KEYS = {
    'time': 'time',
    'flux': 'flux',
    'flux_err': 'flux_err',
}


def parse_filename(filename):
    """Parse TIC ID and sector number from a light curve filename.

    Expected filename format: TIC_<TICID>_<SECTOR>.<ext>. The sector may have
    leading zeros and an optional 'S' prefix (e.g. '06', '6', 'S01', '106').

    Args:
        filename (str): Filename (or path) to parse

    Returns:
        tuple: (tic_id, sector) where tic_id is int and sector is str without
            leading zeros, or (None, None) if the filename cannot be parsed.
    """
    try:
        parts = os.path.basename(filename).split('_')
        tic_id = int(parts[1])
        sector_part = os.path.splitext(parts[2])[0]
        sector = str(int(sector_part.lstrip('Ss')))
        return tic_id, sector
    except (IndexError, ValueError):
        return None, None


def load_light_curve_txt(file_path):
    """Load light curve data from a .txt file.

    The file must contain whitespace-delimited columns in the order
    TIME FLUX FLUX_ERR [PHASE [ECL_MASK]]. A header line is optional and may
    be either '#'-prefixed or bare (detected by attempting to parse the first
    token as a float).

    Args:
        file_path (str): Path to the .txt file

    Returns:
        tuple: (time, flux, flux_err, phase, ecl_mask) where phase and ecl_mask
            are None if the corresponding columns are not present, and ecl_mask
            is boolean (True = in-eclipse).

    Raises:
        ValueError: If the file has fewer than 3 columns
    """
    skiprows = 0
    with open(file_path) as f:
        first = f.readline().strip()
    if first and not first.startswith('#'):
        try:
            float(first.split()[0])
        except ValueError:
            skiprows = 1  # bare (non-commented) header line

    data = np.atleast_2d(np.loadtxt(file_path, skiprows=skiprows))

    if data.shape[1] < 3:
        raise ValueError(
            f"Expected at least 3 columns (TIME FLUX FLUX_ERR) in {file_path}, "
            f"got {data.shape[1]}"
        )

    time = data[:, 0]
    flux = data[:, 1]
    flux_err = data[:, 2]
    phase = data[:, 3] if data.shape[1] > 3 else None
    ecl_mask = data[:, 4].astype(bool) if data.shape[1] > 4 else None

    return time, flux, flux_err, phase, ecl_mask


def load_light_curve_npz(file_path, npz_keys=None):
    """Load light curve data from a .npz file.

    Args:
        file_path (str): Path to the .npz file
        npz_keys (dict, optional): Mapping of logical key names ('time', 'flux',
            'flux_err') to the actual keys in the npz file. Defaults to
            {'time': 'time', 'flux': 'flux', 'flux_err': 'flux_err'}.

    Returns:
        tuple: (time, flux, flux_err, phase, ecl_mask) where phase and ecl_mask
            are None if not present in the file, and ecl_mask is boolean
            (True = in-eclipse).
    """
    if npz_keys is None:
        npz_keys = DEFAULT_NPZ_KEYS

    data = np.load(file_path)
    time = data[npz_keys['time']]
    flux = data[npz_keys['flux']]
    flux_err = data[npz_keys['flux_err']]
    phase = data['phase'] if 'phase' in data else None
    ecl_mask = data['eclipse_mask'].astype(bool) if 'eclipse_mask' in data else None

    return time, flux, flux_err, phase, ecl_mask


def load_light_curve(file_path, npz_keys=None):
    """Load light curve data from a .txt or .npz file based on its extension.

    Args:
        file_path (str): Path to the light curve file
        npz_keys (dict, optional): Key mapping for .npz files (see
            load_light_curve_npz). Ignored for .txt files.

    Returns:
        tuple: (time, flux, flux_err, phase, ecl_mask) where phase and ecl_mask
            are None if not present, and ecl_mask is boolean (True = in-eclipse).

    Raises:
        ValueError: If the file extension is not .txt or .npz
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.npz':
        return load_light_curve_npz(file_path, npz_keys=npz_keys)
    elif ext == '.txt':
        return load_light_curve_txt(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")
