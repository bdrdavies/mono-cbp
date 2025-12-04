"""Eclipse masking for eclipsing binary light curves."""

import numpy as np
import os
import logging
import matplotlib.pyplot as plt
from ..utils.eclipses import get_eclipse_mask, time_to_phase
from ..utils.data import get_row, process_tebc_catalogue
from ..utils import load_catalogue

logger = logging.getLogger('mono_cbp.eclipse_masking')

# Constants
GAP_THRESHOLD = 0.5
POLYFIT_DEGREE = 2
PLOT_FIGSIZE = (10, 5)
PHASE_XLIM_MIN = -0.05
PHASE_XLIM_MAX = 1.05


class EclipseMasker:
    """Class for masking eclipses in EB light curves.

    This class handles the identification and masking of primary and secondary eclipses in TESS light
    curves based on a user-provided catalogue containing binary ephemerides and eclipse parameters.

    Attributes:
        catalogue (pd.DataFrame): User-provided catalogue containing binary ephemerides and eclipse parameters
        data_dir (str): Directory containing light curve data files
    """

    def __init__(self, catalogue, data_dir='./data', TEBC=False):
        """Initialise EclipseMasker.

        Args:
            catalogue (str or pd.DataFrame): Path to or DataFrame of catalogue with eclipse parameters
                (prim_pos, prim_width, sec_pos, sec_width), ephemerides (period, bjd0), and TIC IDs.
            data_dir (str, optional): Directory containing light curve files. Defaults to './data'.
            TEBC (bool, optional): If True, processes TEBC catalogue format with *_2g and *_pf columns
                and converts to standard eclipse parameter columns. If a DataFrame is passed that already
                has standard columns, TEBC processing is skipped. Defaults to False.

        Raises:
            FileNotFoundError: If data_dir does not exist
            NotADirectoryError: If data_dir is not a directory
        """
        if not os.path.exists(data_dir):
            raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
        if not os.path.isdir(data_dir):
            raise NotADirectoryError(f"Path is not a directory: {data_dir}")

        # Load catalogue
        if isinstance(catalogue, str):
            self.catalogue = load_catalogue(catalogue, TEBC=TEBC)
            logger.info(f"Loaded catalogue from {catalogue}")
        else:
            # If DataFrame is passed directly and has TEBC columns, process it only if not already processed
            if TEBC:
                # Check if catalogue has already been processed (has standard eclipse parameter columns)
                required_cols = {'prim_pos', 'prim_width', 'sec_pos', 'sec_width'}
                if not required_cols.issubset(catalogue.columns):
                    catalogue = process_tebc_catalogue(catalogue)
                else:
                    logger.info("Catalogue already has standard eclipse parameter columns, skipping TEBC processing")
            self.catalogue = catalogue

        self.data_dir = data_dir
        logger.info(f"Initialised EclipseMasker with data directory: {data_dir}")

    def _load_npz(self, file_path):
        """Load light curve data from .npz file.

        Args:
            file_path (str): Path to the .npz file

        Returns:
            tuple: (np.ndarray, np.ndarray, np.ndarray, np.ndarray or None, np.ndarray or None)
                - time: Time values from file
                - flux: Flux values from file
                - flux_err: Flux error values from file
                - phase: Binary phase values for each cadence (None if not present in file)
                - eclipse_mask_existing: Existing eclipse mask (None if not present in file)
        """
        npz_data = np.load(file_path)
        time = npz_data['time']
        flux = npz_data['flux']
        flux_err = npz_data['flux_err']
        phase = npz_data.get('phase', None)
        eclipse_mask_existing = npz_data.get('eclipse_mask', None)
        return time, flux, flux_err, phase, eclipse_mask_existing

    def _load_txt(self, file_path):
        """Load light curve data from .txt file.

        Args:
            file_path (str): Path to the .txt file

        Returns:
            tuple: (np.ndarray, np.ndarray, np.ndarray, np.ndarray or None, np.ndarray or None)
                - time: Time values from file
                - flux: Flux values from file
                - flux_err: Flux error values from file
                - phase: Binary phase values (None if not present in file)
                - eclipse_mask_existing: Existing eclipse mask (None if not present in file)
        """
        data = np.loadtxt(file_path, skiprows=1)
        time = data[:, 0]
        flux = data[:, 1]
        flux_err = data[:, 2]

        # Check if phase column exists (column 3)
        if data.shape[1] > 3:
            phase = data[:, 3]
            # Convert eclipse_mask to boolean (stored as 0/1 integers)
            eclipse_mask_existing = data[:, 4].astype(bool) if data.shape[1] > 4 else None
        else:
            phase = None
            eclipse_mask_existing = None

        return time, flux, flux_err, phase, eclipse_mask_existing

    def _combine_eclipse_masks(self, prim_ecl_mask, sec_ecl_mask, phase):
        """Combine primary and secondary eclipse masks.

        Args:
            prim_ecl_mask (np.ndarray): Primary eclipse mask
            sec_ecl_mask (np.ndarray): Secondary eclipse mask
            phase (np.ndarray): Binary phase values (used to initialise empty mask if needed)

        Returns:
            np.ndarray: Combined eclipse mask
        """

        if np.any(prim_ecl_mask) and np.any(sec_ecl_mask):
            result = np.logical_or(prim_ecl_mask, sec_ecl_mask)
        elif np.any(prim_ecl_mask):
            result = prim_ecl_mask
        elif np.any(sec_ecl_mask):
            result = sec_ecl_mask
        else:
            result = np.zeros_like(phase, dtype=bool)

        return result

    def _parse_filename(self, filename):
        """Parse TIC ID and sector number from filename.

        Expected filename format: TIC_<TICID>_<SECTOR>.<EXT>
        Sector is always 2 digits (e.g., 02 or 10) but returned without leading zeros.

        Args:
            filename (str): Filename to parse

        Returns:
            tuple: (tic_id, sector) where tic_id is int and sector is str without leading zeros

        Raises:
            ValueError: If filename format is invalid
        """
        try:
            parts = filename.split('_')
            tic_id = int(parts[1])
            sector_part = parts[2].split('.')[0]
            sector_num = int(sector_part)  # Convert to int to remove leading zero
            sector = str(sector_num)
            return tic_id, sector
        except (ValueError, IndexError) as e:
            raise ValueError(f"Failed to parse TIC ID and sector from filename '{filename}': {e}")

    def _detrend_out_of_eclipse_data(self, data):
        """Detrend out-of-eclipse flux using segmented polynomial fitting.

        Applies polynomial detrending only to out-of-eclipse (OOE) data. Data is split
        into segments by time gaps, and a quadratic polynomial is fit independently to
        OOE points within each segment. Then flux is divided by the fitted trend.
        This preserves eclipse features while removing out-of-eclipse variations (for plotting purposes).

        Args:
            data (np.ndarray): Data array with shape (N, 5) and columns:
                [time, flux, flux_err, phase, eclipse_mask]
                eclipse_mask: boolean where True = in-eclipse

        Returns:
            np.ndarray: Modified data array with detrended flux (column 1).
                        Returns data unchanged if <5 columns or no OOE data to fit.
        """

        if data.shape[1] != 5 or not np.any(data[:, -1] == False):
            return data

        all_time = data[:, 0]
        all_flux = data[:, 1]
        # Convert eclipse mask to boolean
        ecl_mask_bool = data[:, -1].astype(bool)
        ooe_time = data[~ecl_mask_bool, 0]
        all_gaps = np.where(np.diff(all_time) > GAP_THRESHOLD)[0]
        all_segments = np.split(np.arange(len(all_time)), all_gaps + 1)

        for all_seg in all_segments:
            seg_time = all_time[all_seg]
            ooe_mask = np.isin(seg_time, ooe_time)
            ooe_seg_time = seg_time[ooe_mask]
            ooe_seg_flux = all_flux[all_seg][ooe_mask]

            if len(ooe_seg_time) < 2:
                continue

            # Check for NaN values in out-of-eclipse flux
            if np.any(np.isnan(ooe_seg_flux)):
                continue

            p = np.polyfit(ooe_seg_time, ooe_seg_flux, POLYFIT_DEGREE)

            # Check if polyfit returned NaN
            if np.any(np.isnan(p)):
                continue

            model = np.polyval(p, seg_time)
            model[model == 0] = 1.0
            data[all_seg, 1] /= model

        return data

    def mask_all(self, force=False):
        """Mask eclipses in all .txt and .npz files in the data directory.

        Skips files that cannot be processed (e.g., missing TIC ID in catalogue, unsupported format).

        Args:
            force (bool, optional): If True, recalculate mask even if it already exists. Defaults to False.
        """
        files = [f for f in os.listdir(self.data_dir) if f.endswith('.txt') or f.endswith('.npz')]
        logger.info(f"Processing {len(files)} files in {self.data_dir}")

        for file in files:
            try:
                self.mask_file(file, force=force)
            except (FileNotFoundError, ValueError) as e:
                logger.warning(f"Skipping {file}: {e}")

    def mask_file(self, file, force=False):
        """Mask eclipses in a single light curve file (.txt or .npz).

        Loads light curve data (time, flux, flux_err) and optionally phase. If phase is not provided,
        it is calculated from the catalogue ephemeris. Computes eclipse masks for primary and secondary
        eclipses and saves the results with an appended eclipse_mask column.

        Args:
            file (str): Filename (not full path) to process. Filename should be in format 'TIC_<TICID>_<sector>.<ext>'
            force (bool, optional): If True, recalculate mask even if it already exists. Defaults to False.

        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If file format is unsupported or TIC ID not found in catalogue
        """
        file_path = os.path.join(self.data_dir, file)
        file_ext = os.path.splitext(file)[1]

        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        # Load data based on file format
        if file_ext == '.npz':
            time, flux, flux_err, phase, eclipse_mask_existing = self._load_npz(file_path)
        elif file_ext == '.txt':
            time, flux, flux_err, phase, eclipse_mask_existing = self._load_txt(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file}")

        # Extract TIC ID and sector from filename
        tic_id, _ = self._parse_filename(file)

        # Get eclipse parameters from catalogue
        row = get_row(self.catalogue, tic_id)
        if row is None:
            raise ValueError(f"TIC {tic_id} not found in catalogue.")

        # Calculate phase from ephemeris if not provided in input file
        if phase is None:
            phase = time_to_phase(time, row['period'], row['bjd0'])

        # Check if mask already exists
        if eclipse_mask_existing is not None and not force:
            logger.info(f"Eclipse mask already exists for {file}")
            return

        # Create eclipse masks
        prim_ecl_mask = get_eclipse_mask(phase, row['prim_pos'], row['prim_width'])
        sec_ecl_mask = get_eclipse_mask(phase, row['sec_pos'], row['sec_width'])

        # Combine masks
        eclipse_mask = self._combine_eclipse_masks(prim_ecl_mask, sec_ecl_mask, phase)

        # Save based on format
        if file_ext == '.npz':
            # Save as npz
            np.savez(file_path,
                     time=time,
                     flux=flux,
                     flux_err=flux_err,
                     phase=phase,
                     eclipse_mask=eclipse_mask)
            logger.info(f"Saved eclipse mask to {file}")
        else:
            # Save as txt
            data = np.column_stack((time, flux, flux_err, phase, eclipse_mask.astype(int)))
            np.savetxt(file_path, data, header="TIME FLUX FLUX_ERR PHASE ECL_MASK")
            logger.info(f"Saved eclipse mask to {file}")

    def plot_bin_phase_fold(self, tic_id, save_fig=False, save_path='.'):
        """Plot phase-folded light curve with in-eclipse data highlighted.

        Loads and combines all data files for the given TIC ID. If phase is not provided in the input files,
        it is calculated from the catalogue ephemeris. Plots in-eclipse and out-of-eclipse data with different colors.
        If eclipse mask exists, applies polynomial detrending to out-of-eclipse data before plotting.

        Args:
            tic_id (int): TIC ID to plot
            save_fig (bool, optional): Whether to save the figure. Defaults to False.
            save_path (str, optional): Directory to save figure. Defaults to working directory ('.').

        Raises:
            FileNotFoundError: If no files found for the given TIC ID
            ValueError: If no valid data can be loaded for plotting
        """
        # Find all files for this TIC (both .txt and .npz)
        files = [f for f in os.listdir(self.data_dir)
                 if (f.endswith('.txt') or f.endswith('.npz')) and f.split('_')[1] == str(tic_id)]

        if not files:
            raise FileNotFoundError(f"No files found for TIC {tic_id}")

        # Extract sectors from filenames
        sectors = [self._parse_filename(f)[1] for f in files]
        sectors_str = ','.join(sorted(sectors, key=int))

        # Get period from catalogue
        row = get_row(self.catalogue, tic_id)
        period = row['period'] if row is not None else 'unknown'

        # Load and combine all data
        all_data = []
        for file in files:
            file_path = os.path.join(self.data_dir, file)

            # Load data based on file format
            if file.endswith('.npz'):
                time, flux, flux_err, phase, eclipse_mask = self._load_npz(file_path)
            else:
                time, flux, flux_err, phase, eclipse_mask = self._load_txt(file_path)

            # Calculate phase from ephemeris if not provided
            if phase is None:
                phase = time_to_phase(time, row['period'], row['bjd0'])

            # Calculate eclipse mask if not provided
            if eclipse_mask is None:
                prim_ecl_mask = get_eclipse_mask(phase, row['prim_pos'], row['prim_width'])
                sec_ecl_mask = get_eclipse_mask(phase, row['sec_pos'], row['sec_width'])
                eclipse_mask = self._combine_eclipse_masks(prim_ecl_mask, sec_ecl_mask, phase)
            else:
                # Ensure eclipse_mask is 1D boolean array
                eclipse_mask = np.asarray(eclipse_mask).ravel().astype(bool)

            # Create array in same format: TIME FLUX FLUX_ERR PHASE ECL_MASK
            data = np.column_stack((time, flux, flux_err, phase, eclipse_mask))
            all_data.append(data)

        if not all_data:
            raise ValueError(f"No valid data to plot for TIC {tic_id}")

        data = np.vstack(all_data)

        data = data[np.argsort(data[:, 0])]

        # Detrend out-of-eclipse data if mask exists
        data = self._detrend_out_of_eclipse_data(data)

        # Create plot
        fig = plt.figure(figsize=PLOT_FIGSIZE)

        if data.shape[1] == 5:
            # Ensure eclipse mask is boolean
            eclipse_mask = data[:, -1].astype(bool)

            # Plot with eclipse mask distinction
            plt.errorbar(data[eclipse_mask, 3], data[eclipse_mask, 1],
                        yerr=data[eclipse_mask, 2], fmt='.', color='red',
                        alpha=0.5, label='In-eclipse')
            plt.errorbar(data[~eclipse_mask, 3], data[~eclipse_mask, 1],
                        yerr=data[~eclipse_mask, 2], fmt='.', color='navy',
                        alpha=0.5, label='Out-of-eclipse')
        else:
            # Plot without eclipse mask distinction
            plt.errorbar(data[:, 3], data[:, 1], yerr=data[:, 2], fmt='.', color='black', alpha=0.5)

        plt.xlabel('Phase')
        plt.ylabel('Normalized Flux')
        period_str = f'{period:.2f} days' if isinstance(period, (int, float)) else period
        plt.title(f'TIC {tic_id} (Sectors: {sectors_str})\nPeriod: {period_str}')
        # Set x-axis limits to ensure all phase values are visible, including those at boundaries
        plt.xlim(PHASE_XLIM_MIN, PHASE_XLIM_MAX)
        plt.legend()

        if save_fig:
            output_path = f"{save_path}/TIC_{tic_id}_phase-fold.png"
            plt.savefig(output_path)
            logger.info(f"Saved phase-fold plot to {output_path}")
            plt.close(fig)
        else:
            plt.show()