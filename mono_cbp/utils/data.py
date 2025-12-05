"""Collection of utility functions for loading and pre-processing TESS data."""

import numpy as np
import pandas as pd
import os
import logging
import lightkurve as lk
import astropy.units as u
from .eclipses import time_to_phase

logger = logging.getLogger('mono_cbp.utils.data')


def _is_valid(val, max_width=None):
    """Check if a value is valid (not null and not zero).

    Args:
        val: Value to check
        max_width (float, optional): Maximum allowed width value (for eclipse widths)

    Returns:
        bool: True if valid, False otherwise
    """
    if not pd.notnull(val) or val == 0:
        return False

    # Apply max width check if provided
    if max_width is not None and val > max_width:
        return False

    return True


def choose_eclipse_params(row, max_eclipse_width=0.2):
    """Select eclipse parameters using only 2g (Double Gaussian) values.

    This function only uses 2g values for eclipse parameter selection, ignoring
    all polyfit (pf) values.

    Args:
        row (pd.Series): Row from the catalogue
        max_eclipse_width (float, optional): Maximum allowed eclipse width in phase units. Default 0.2.

    Returns:
        dict: Dictionary with keys 'prim_pos', 'prim_width', 'sec_pos', 'sec_width'
    """
    result = {
        'prim_pos': 0,
        'prim_width': 0,
        'sec_pos': 0,
        'sec_width': 0
    }

    # Primary eclipse (2g only)
    prim_2g_valid = (_is_valid(row.get('prim_pos_2g')) and
                     _is_valid(row.get('prim_width_2g'), max_width=max_eclipse_width))

    if prim_2g_valid:
        result['prim_pos'] = row['prim_pos_2g']
        result['prim_width'] = row['prim_width_2g']
    elif _is_valid(row.get('prim_pos_2g')):
        result['prim_pos'] = row['prim_pos_2g']

    if _is_valid(row.get('prim_width_2g'), max_width=max_eclipse_width):
        result['prim_width'] = row['prim_width_2g']

    # Secondary eclipse (2g only)
    sec_2g_valid = (_is_valid(row.get('sec_pos_2g')) and
                    _is_valid(row.get('sec_width_2g'), max_width=max_eclipse_width))

    if sec_2g_valid:
        result['sec_pos'] = row['sec_pos_2g']
        result['sec_width'] = row['sec_width_2g']
    elif _is_valid(row.get('sec_pos_2g')):
        result['sec_pos'] = row['sec_pos_2g']

    if _is_valid(row.get('sec_width_2g'), max_width=max_eclipse_width):
        result['sec_width'] = row['sec_width_2g']

    return result


def choose_pf_2g(row, col_2g, col_pf):
    """Helper function to choose between double Gaussian and polyfit values.

    Only applicable to TEBC.

    Args:
        row (pd.Series): Row from the catalogue
        col_2g (str): Column name for double Gaussian value
        col_pf (str): Column name for polyfit value

    Returns:
        float: Selected value
    """
    val_2g = row[col_2g]
    return val_2g if pd.notnull(val_2g) and val_2g != 0 else row[col_pf]


def process_tebc_catalogue(catalogue):
    """Process a TEBC-format catalogue DataFrame to extract eclipse parameters.

    Converts TEBC-format columns (*_2g and *_pf variants) to standard eclipse parameters
    (prim_pos, prim_width, sec_pos, sec_width) using the choose_eclipse_params function.

    Args:
        catalogue (pd.DataFrame): DataFrame with TEBC columns

    Returns:
        pd.DataFrame: Modified catalogue with standard eclipse parameter columns added

    Raises:
        ValueError: If catalogue is not a DataFrame
    """
    if not isinstance(catalogue, pd.DataFrame):
        raise ValueError("Catalogue must be a DataFrame")

    catalogue = catalogue.copy()
    eclipse_params = catalogue.apply(choose_eclipse_params, axis=1)
    catalogue['prim_pos'] = eclipse_params.apply(lambda x: x['prim_pos'])
    catalogue['prim_width'] = eclipse_params.apply(lambda x: x['prim_width'])
    catalogue['sec_pos'] = eclipse_params.apply(lambda x: x['sec_pos'])
    catalogue['sec_width'] = eclipse_params.apply(lambda x: x['sec_width'])
    logger.info("Processed TEBC catalogue: eclipse parameters selected from 2g columns")
    return catalogue


def load_catalogue(path, TEBC=False):
    """Load a catalogue of TESS eclipsing binaries from a CSV or TXT file.

    The catalogue must contain the following columns:
        - tess_id: TIC IDs of the targets
        - period: Orbital period of the targets
        - bjd0: Reference epoch for binary ephemerides
        - sectors: TESS sectors that the targets were observed in
        - prim_pos: Primary eclipse position (orbital phase)
        - prim_width: Primary eclipse width (orbital phase)
        - sec_pos: Secondary eclipse position (orbital phase)
        - sec_width: Secondary eclipse width (orbital phase)

    Args:
        path (str): Path to the catalogue file (.csv or .txt).
        TEBC (bool, optional): If True, use TEBC double Gaussian/polyfit logic. Requires
            additional columns: prim_pos_2g, prim_pos_pf, prim_width_2g, prim_width_pf,
            sec_pos_2g, sec_pos_pf, sec_width_2g, sec_width_pf. Default: False.

    Returns:
        pd.DataFrame: Catalogue with required columns.

    Raises:
        TypeError: If path is not a string
        FileNotFoundError: If catalogue file doesn't exist
        ValueError: If file format is unsupported or required columns are missing
    """
    if not isinstance(path, str):
        raise TypeError("Path must be a string.")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Catalogue file not found at {path}")

    if path.endswith('.txt'):
        catalogue = pd.read_fwf(path)
    elif path.endswith('.csv'):
        catalogue = pd.read_csv(path)
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .txt file.")

    required_columns = [
        'tess_id', 'period', 'bjd0', 'sectors',
        'prim_pos', 'prim_width', 'sec_pos', 'sec_width'
    ]

    # If TEBC format, create eclipse width/position columns from 2g/pf values
    if TEBC:
        eclipse_params = catalogue.apply(choose_eclipse_params, axis=1)
        catalogue['prim_pos'] = eclipse_params.apply(lambda x: x['prim_pos'])
        catalogue['prim_width'] = eclipse_params.apply(lambda x: x['prim_width'])
        catalogue['sec_pos'] = eclipse_params.apply(lambda x: x['sec_pos'])
        catalogue['sec_width'] = eclipse_params.apply(lambda x: x['sec_width'])

    # Check for required columns (after creating them if TEBC=True)
    for col in required_columns:
        if col not in catalogue.columns:
            raise ValueError(f"Catalogue is missing required column: {col}")

    logger.info(f"Loaded catalogue with {len(catalogue)} targets from {path}")
    return catalogue[required_columns]


def get_row(catalogue, tic_id):
    """Get the row in a catalogue DataFrame corresponding to a specific TIC ID.

    Args:
        catalogue (pd.DataFrame or pd.Series): Catalogue data
        tic_id (int): TIC ID to search for

    Returns:
        pd.Series or None: Row data for the TIC ID, or None if not found
    """
    if isinstance(catalogue, pd.Series):
        return catalogue if catalogue['tess_id'] == tic_id else None
    else:
        row = catalogue[catalogue['tess_id'] == tic_id]
        return row.iloc[0] if not row.empty else None


def bin_to_long_cadence(time, flux, flux_err):
    """Bins faster cadence data to 30 minute cadence.

    Args:
        time (array): Time values of the light curve to bin
        flux (array): Flux values of the light curve to bin
        flux_err (array): Flux error values of the light curve to bin

    Returns:
        tuple: time, flux, and flux_err binned to 30 minutes
    """
    # Create a lightkurve LightCurve object
    lc = lk.LightCurve(time=np.nan_to_num(time), flux=np.nan_to_num(flux), flux_err=np.nan_to_num(flux_err))

    # Set NaN values
    lc.time.value[lc.time.value == 0] = np.nan
    lc.flux.value[lc.flux.value == 0] = np.nan
    lc.flux_err.value[lc.flux_err.value == 0] = np.nan

    # Bin the light curve to 30 minute cadence
    lc_binned = lc.bin(time_bin_size=30*u.minute)

    # Extract binned values
    time = lc_binned.time.value
    flux = lc_binned.flux.value
    flux_err = lc_binned.flux_err.value

    return time, flux, flux_err


def _to_plain(arr):
    """Convert MaskedNDArray to plain ndarray with NaNs for masked values.

    Args:
        arr: Array to convert

    Returns:
        np.ndarray: Plain numpy array
    """
    return arr.filled(np.nan) if hasattr(arr, "filled") else np.asarray(arr)


def lc_to_txt(cat, lc, output_path='../data'):
    """Saves a lightkurve LightCurve object to a text file in the standard mono-cbp format.

    Args:
        cat (pd.DataFrame or pd.Series): Catalogue data containing period and bjd0
        lc (lightkurve.LightCurve): LightCurve object to save
        output_path (str, optional): Directory to save the file. Defaults to './data'.
    """
    # Create filename
    tic_id = lc.meta.get("TARGETID", "unknown")
    sector = lc.meta.get("SECTOR", "unknown")
    filename = f"TIC_{tic_id}_{sector:02d}.txt"
    output_file = os.path.join(output_path, filename)

    # If output file already exists, skip saving
    if os.path.exists(output_file):
        logger.debug(f"File {output_file} already exists, skipping")
        return

    time = lc.time.value
    raw_flux = lc.flux.value
    raw_flux_err = lc.flux_err.value

    # Obtain CROWDSAP and FLFRCSAP values for contamination correction
    CROWDSAP = lc.meta.get("CROWDSAP")
    FLFRCSAP = lc.meta.get("FLFRCSAP")

    # Set flux values of 0 to NaN
    raw_flux[raw_flux == 0] = np.nan
    raw_flux_err[raw_flux == 0] = np.nan

    # Perform the FF and CM corrections
    median_flux = np.nanmedian(raw_flux)
    excess_flux = (1 - CROWDSAP) * median_flux
    flux_removed = raw_flux - excess_flux
    flux = flux_removed / FLFRCSAP
    flux_err = raw_flux_err / FLFRCSAP

    # Normalise flux and flux error values
    flux_err /= np.nanmedian(flux)
    flux /= np.nanmedian(flux)

    # Calculate phase values
    if isinstance(cat, pd.Series):
        period = cat['period']
        bjd0 = cat['bjd0']
    else:
        period = cat[cat['tess_id'] == tic_id]['period'].values[0]
        bjd0 = cat[cat['tess_id'] == tic_id]['bjd0'].values[0]
    phase = time_to_phase(time, period=period, t0=bjd0)

    # Bin to 30 minute cadence if data is at 2 minute cadence
    if np.nanmedian(np.diff(time)) < 0.0208:  # 0.0208 days = 30 minutes
        time, flux, flux_err = bin_to_long_cadence(time, flux, flux_err)
        phase = time_to_phase(time, period=period, t0=bjd0)

    # Filter NaN values
    nan_mask = ~np.isnan(time * flux * flux_err * phase)
    time = time[nan_mask]
    flux = flux[nan_mask]
    flux_err = flux_err[nan_mask]
    phase = phase[nan_mask]

    # Save the extracted data to a text file
    np.savetxt(
        output_file,
        np.column_stack([
            _to_plain(time.astype(float)),
            _to_plain(flux.astype(float)),
            _to_plain(flux_err.astype(float)),
            _to_plain(phase.astype(float))
        ]),
        header="TIME FLUX FLUX_ERR PHASE ECL_MASK",
    )
    logger.info(f"Saved light curve data to {output_file}")


def catalogue_to_lc_files(cat, output_path='../data'):
    """Creates light curve files for each TIC ID in the catalogue in the standard mono-cbp format.

    Args:
        cat (pd.DataFrame or pd.Series): The catalogue data
        output_path (str, optional): The directory to save the downloaded lc.fits files. Defaults to './data'.
    """
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        logger.info(f"Created output directory: {output_path}")

    if isinstance(cat, pd.Series):
        tic_id = cat['tess_id']
        sectors = str(cat['sectors']).split(',')
        for sector in sectors:
            sector = int(sector)
            # Check if a file already exists for this TIC ID and sector in the output directory
            existing_files = [f for f in os.listdir(output_path) if f.endswith('.txt') and f.split('_')[1] == str(tic_id) and f.split('_')[2].replace('.txt', '') == f"{sector:02d}"]
            if existing_files:
                logger.info(f"TIC_{tic_id}_{sector:02d}.txt already exists in {output_path}. Skipping download.")
                continue
            try:
                search = lk.search_lightcurve(f"TIC {tic_id}", mission="TESS", sector=sector, author="TESS-SPOC")
                if len(search) == 0:
                    logger.warning(f"No TESS-SPOC light curve found for TIC {tic_id} in sector {sector}.")
                    continue
                lc = search.download(flux_column='sap_flux', quality_bitmask='hard')
                if lc is None:
                    logger.warning(f"Download failed from MAST for TIC {tic_id} in sector {sector}.")
                    continue
                lc_to_txt(cat, lc, output_path=output_path)
            except Exception as e:
                logger.error(f"Failed to download data for TIC {tic_id} in sector {sector}: {e}")

    else:
        for i in range(len(cat)):
            tic_id = cat['tess_id'].iloc[i]
            sectors = str(cat['sectors'].iloc[i]).split(',')
            for sector in sectors:
                sector = int(sector)
                # Check if a file already exists for this TIC ID and sector in the output directory
                existing_files = [f for f in os.listdir(output_path) if f.endswith('.txt') and f.split('_')[1] == str(tic_id) and f.split('_')[2].replace('.txt', '') == f"{sector:02d}"]
                if existing_files:
                    logger.info(f"TIC_{tic_id}_{sector:02d}.txt already exists in {output_path}. Skipping download.")
                    continue
                try:
                    search = lk.search_lightcurve(f"TIC {tic_id}", mission="TESS", sector=sector, author="TESS-SPOC")
                    if len(search) == 0:
                        logger.warning(f"No TESS-SPOC light curve found for TIC {tic_id} in sector {sector}.")
                        continue
                    lc = search.download(flux_column='sap_flux', quality_bitmask='hard')
                    if lc is None:
                        logger.warning(f"Download failed from MAST for TIC {tic_id} in sector {sector}.")
                        continue
                    lc_to_txt(cat, lc, output_path=output_path)
                except Exception as e:
                    logger.error(f"Failed to download data for TIC {tic_id} in sector {sector}: {e}")
